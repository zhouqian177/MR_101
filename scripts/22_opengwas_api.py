#!/usr/bin/env python3
# 22_opengwas_api.py — OpenGWAS API 客户端（Python 通道）
#
# 背景: R 的 libcurl 7.61.1 过旧无法通过代理 TLS 握手(SSL_ERROR_SYSCALL),
#       Python requests 可正常访问 api.opengwas.io。故用 Python 完成 API 数据获取,
#       结果存为 CSV/JSON, 供 R 分析层读取。
#
# 用法:
#   python3 scripts/22_opengwas_api.py gwasinfo --ids ieu-b-110,ieu-a-7
#   python3 scripts/22_opengwas_api.py associations --id ieu-a-7 --snps rs1800588,rs599839
#   python3 scripts/22_opengwas_api.py tophits --id ieu-b-110 --p1 5e-8 --out 02.analysis/opengwas/online_instruments.csv
#   python3 scripts/22_opengwas_api.py search --trait "LDL" --out 02.analysis/opengwas/dataset_search.csv
#
# token: 从 .env 读取 OpenGWAS_JWT（不打印）

import argparse, json, os, sys, time
import requests

API = "https://api.opengwas.io/api/"

def get_token():
    # 优先环境变量, 否则读取 .env（脚本项目根目录）
    tok = os.environ.get("OPENGWAS_JWT", "")
    if not tok:
        env_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".env")
        if os.path.exists(env_path):
            with open(env_path) as f:
                for line in f:
                    line = line.strip().strip('\r')
                    if line.startswith("OpenGWAS_JWT="):
                        tok = line.split("=", 1)[1].strip().strip('"')
    if not tok:
        sys.exit("ERROR: 未找到 OPENGWAS_JWT token（设置环境变量或 .env 文件）")
    return tok

def headers():
    return {"Authorization": "Bearer " + get_token()}

def api_get(path, params=None, retries=3):
    for i in range(retries):
        try:
            r = requests.get(API + path, params=params, headers=headers(), timeout=120)
            if r.status_code == 200:
                return r.json()
            if r.status_code == 429:
                wait = int(r.headers.get("Retry-After", 30))
                print(f"  429 限流, 等待 {wait}s ...", file=sys.stderr)
                time.sleep(wait); continue
            sys.exit(f"ERROR: {path} -> HTTP {r.status_code}: {r.text[:200]}")
        except requests.exceptions.RequestException as e:
            print(f"  重试 {i+1}: {e}", file=sys.stderr); time.sleep(2)
    sys.exit("ERROR: API 请求多次失败")

def api_post(path, payload, retries=3):
    for i in range(retries):
        try:
            r = requests.post(API + path, json=payload, headers=headers(), timeout=300)
            if r.status_code == 200:
                return r.json()
            if r.status_code == 429:
                wait = int(r.headers.get("Retry-After", 30))
                print(f"  429 限流, 等待 {wait}s ...", file=sys.stderr)
                time.sleep(wait); continue
            sys.exit(f"ERROR: POST {path} -> HTTP {r.status_code}: {r.text[:200]}")
        except requests.exceptions.RequestException as e:
            print(f"  重试 {i+1}: {e}", file=sys.stderr); time.sleep(2)
    sys.exit("ERROR: API 请求多次失败")

def cmd_gwasinfo(args):
    ids = args.ids.split(",")
    # 注意: gwasinfo 端点忽略 id 参数, 返回全量数据集 dict(id->info), 本地过滤
    d = api_get("gwasinfo", {})
    for i in ids:
        x = d.get(i)
        if not x:
            print(f"{i}\t未找到"); continue
        print(f"{x['id']}\t{x.get('trait','')}\t"
              f"ncase={x.get('ncase','')}\tncontrol={x.get('ncontrol','')}\tnsnp={x.get('nsnp','')}\t"
              f"population={x.get('population','')}\tbuild={x.get('buildnumber','')}")

def cmd_search(args):
    # /api/gwasinfo 支持 trait 过滤; 这里遍历常见字段做关键词检索
    d = api_get("gwasinfo", {})
    hits = [x for x in d if args.trait.lower() in x.get("trait", "").lower() or
            args.trait.lower() in x.get("id", "").lower()]
    print(f"匹配 {len(hits)} 个数据集")
    out = args.out or "02.analysis/opengwas/dataset_search.csv"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as f:
        f.write("id\ttrait\tncase\tncontrol\tnsnp\tpopulation\n")
        for x in hits:
            f.write(f"{x['id']}\t{x.get('trait','')}\t{x.get('ncase','')}\t"
                    f"{x.get('ncontrol','')}\t{x.get('nsnp','')}\t{x.get('population','')}\n")
    print(f"已保存 {out}")

def cmd_associations(args):
    snps = args.snps.split(",")
    # API 限制: N(id) * N(variant) <= 64, 超过则自动分批查询
    BATCH = 60
    all_rows = []
    for i in range(0, len(snps), BATCH):
        chunk = snps[i:i + BATCH]
        print(f"  批次 {i // BATCH + 1}: {len(chunk)} 个 SNP ...", file=sys.stderr)
        d = api_post("associations", {"variant": chunk, "id": args.id})
        all_rows.extend(d)
    print(f"查询 {len(snps)} 个 SNP × {args.id}: 返回 {len(all_rows)} 条")
    if args.out:
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        with open(args.out, "w") as f:
            f.write("rsid\tea\tnea\teaf\tbeta\tse\tp\n")
            for a in all_rows:
                f.write(f"{a.get('rsid','')}\t{a.get('ea','')}\t{a.get('nea','')}\t"
                        f"{a.get('eaf','')}\t{a.get('beta','')}\t{a.get('se','')}\t{a.get('p','')}\n")
        print(f"已保存 {args.out}")
    for a in all_rows:
        print(f"{a.get('rsid','')}\t{a.get('ea','')}/{a.get('nea','')}\t"
              f"beta={a.get('beta','')}\tse={a.get('se','')}\tp={a.get('p','')}\t"
              f"eaf={a.get('eaf','')}")

def cmd_tophits(args):
    # 注意: tophits 端点为 POST（GET 返回 405）
    d = api_post("tophits", {"id": args.id, "pval": float(args.p1)})
    out = args.out or f"02.analysis/opengwas/online_{args.id}_instruments.csv"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as f:
        f.write("rsid\tchr\tpos\tea\tnea\teaf\tbeta\tse\tp\n")
        for x in d:
            f.write(f"{x.get('rsid','')}\t{x.get('chr','')}\t{x.get('position','')}\t"
                    f"{x.get('ea','')}\t{x.get('nea','')}\t{x.get('eaf','')}\t"
                    f"{x.get('beta','')}\t{x.get('se','')}\t{x.get('p','')}\n")
    print(f"tophits(P<{args.p1}): {len(d)} 个 SNP, 已保存 {out}")

def cmd_phewas(args):
    # PheWAS: 单个变异跨全部数据集的关联扫描 (POST /api/phewas)
    d = api_post("phewas", {"variant": args.snps, "pval": float(args.p1)})
    out = args.out or "02.analysis/opengwas/phewas.csv"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as f:
        f.write("id\ttrait\tchr\tpos\trsid\tea\tnea\teaf\tbeta\tse\tp\tn\n")
        for a in d:
            f.write(f"{a.get('id','')}\t{a.get('trait','')}\t{a.get('chr','')}\t"
                    f"{a.get('position','')}\t{a.get('rsid','')}\t{a.get('ea','')}\t"
                    f"{a.get('nea','')}\t{a.get('eaf','')}\t{a.get('beta','')}\t"
                    f"{a.get('se','')}\t{a.get('p','')}\t{a.get('n','')}\n")
    print(f"PheWAS({args.snps}, P<{args.p1}): {len(d)} 条关联, 已保存 {out}")

def main():
    p = argparse.ArgumentParser(description="OpenGWAS API 客户端")
    sub = p.add_subparsers(dest="cmd", required=True)
    g1 = sub.add_parser("gwasinfo"); g1.add_argument("--ids", required=True); g1.set_defaults(fn=cmd_gwasinfo)
    g2 = sub.add_parser("search"); g2.add_argument("--trait", required=True); g2.add_argument("--out", default=None); g2.set_defaults(fn=cmd_search)
    g3 = sub.add_parser("associations"); g3.add_argument("--id", required=True); g3.add_argument("--snps", required=True); g3.add_argument("--out", default=None); g3.set_defaults(fn=cmd_associations)
    g4 = sub.add_parser("tophits"); g4.add_argument("--id", required=True); g4.add_argument("--p1", default="5e-8"); g4.add_argument("--out", default=None); g4.set_defaults(fn=cmd_tophits)
    g5 = sub.add_parser("phewas"); g5.add_argument("--snps", required=True); g5.add_argument("--p1", default="5e-8"); g5.add_argument("--out", default=None); g5.set_defaults(fn=cmd_phewas)
    args = p.parse_args()
    args.fn(args)

if __name__ == "__main__":
    main()
