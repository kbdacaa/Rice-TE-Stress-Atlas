"""自动下载GEO表达数据集
- 从stress_expression_catalog.json读取416个数据集
- 对含GEO编号的(GSE*)，通过NCBI edirect自动下载series matrix
- 对SRA数据(PRJNA/SRP)，输出wget命令列表供手动下载
- 仅下载补充文件(processed data)，不下载原始FASTQ
"""
import json, os, subprocess, time, sys

CATALOG = r"D:\project\AutoResearch-rice-T2T\data\stress_expression_catalog.json"
OUT_DIR = r"D:\fanjiyinzu\NCBI_expression"

def get_geo_suppl(geo_id):
    """Download supplementary files from GEO using NCBI edirect"""
    # GEO series matrix URL pattern
    # https://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnnnn/GSE76415/matrix/
    series_dir = f"GSE{geo_id[-3:]}nnn"
    # Actually use the standard pattern
    gse_num = geo_id.replace('GSE','')

    # Try to get the series matrix via NCBI FTP
    ftp_base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{geo_id[:5]}nnn/{geo_id}"

    urls = [
        f"{ftp_base}/matrix/{geo_id}_series_matrix.txt.gz",
        f"{ftp_base}/suppl/{geo_id}_FPKM.txt.gz",
        f"{ftp_base}/suppl/{geo_id}_fpkm.txt.gz",
        f"{ftp_base}/suppl/{geo_id}_counts.txt.gz",
        f"{ftp_base}/suppl/{geo_id}_normalized.txt.gz",
        f"{ftp_base}/suppl/{geo_id}_expr.txt.gz",
    ]
    return urls

def try_download(url, out_path):
    """Try to download a URL"""
    try:
        result = subprocess.run(
            ['curl', '-s', '-L', '--connect-timeout', '10', '--max-time', '60',
             '-o', out_path, url],
            capture_output=True, text=True, timeout=70
        )
        if os.path.exists(out_path) and os.path.getsize(out_path) > 500:
            return True
        if os.path.exists(out_path):
            os.remove(out_path)
        return False
    except Exception as e:
        return False

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    with open(CATALOG, 'r', encoding='utf-8') as f:
        cat = json.load(f)

    datasets = cat['datasets']

    # Filter by download type
    geo_datasets = [d for d in datasets if d.get('geo') and d['geo'][0] not in
                    ['GSE76415','GSE182563','GSE92989','GSE41647','GSE60287','GSE98455','GSE79011']]
    sra_datasets = [d for d in datasets if (not d.get('geo') or not d['geo']) and d.get('sra')]

    print(f"待下载GEO: {len(geo_datasets)}个")
    print(f"待下载SRA(仅输出命令): {len(sra_datasets)}个")
    print(f"输出目录: {OUT_DIR}")
    print()

    # ─── Download GEO processed data ───
    success = 0
    failed = 0
    skipped = 0

    for i, ds in enumerate(geo_datasets):  # ALL GEO datasets
        geo_id = ds['geo'][0]
        title = ds.get('title', '')[:60]
        stresses = '+'.join(ds.get('stresses', []))

        out_subdir = os.path.join(OUT_DIR, geo_id)

        # Check if already downloaded
        if os.path.exists(out_subdir) and os.listdir(out_subdir):
            skipped += 1
            continue

        os.makedirs(out_subdir, exist_ok=True)

        # Try to download series matrix (universal format)
        ftp_base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{geo_id[:5]}nnn/{geo_id}"
        matrix_url = f"{ftp_base}/matrix/{geo_id}_series_matrix.txt.gz"
        out_matrix = os.path.join(out_subdir, f"{geo_id}_series_matrix.txt.gz")

        print(f"[{i+1}/{len(geo_datasets[:50])}] {geo_id} [{stresses}] {title}")

        if try_download(matrix_url, out_matrix):
            size_kb = os.path.getsize(out_matrix) / 1024
            print(f"  OK series_matrix ({size_kb:.0f}KB)")
            success += 1
        else:
            # Try suppl directory
            suppl_url = f"{ftp_base}/suppl/"
            # Download filelist
            filelist_path = os.path.join(out_subdir, "filelist.txt")
            if try_download(f"{ftp_base}/suppl/filelist.txt", filelist_path):
                print(f"  Got filelist, checking for processed data...")
                with open(filelist_path, 'r') as f:
                    for line in f:
                        if any(x in line.lower() for x in ['fpkm','rpkm','count','norm','expr','tpm']):
                            # Extract filename
                            parts = line.strip().split('\t')
                            if parts:
                                fname = parts[0].strip()
                                if fname.endswith('.gz'):
                                    suppl_furl = f"{ftp_base}/suppl/{fname}"
                                    suppl_out = os.path.join(out_subdir, fname)
                                    print(f"    Trying: {fname}")
                                    if try_download(suppl_furl, suppl_out):
                                        print(f"    OK {fname}")
                            break
                success += 1
            else:
                print(f"  FAILED (no matrix, no filelist)")
                failed += 1

        # Be nice to NCBI
        time.sleep(0.5)

    print(f"\n━━━ GEO下载完成 ━━━")
    print(f"成功: {success}, 失败: {failed}, 已跳过: {skipped}")
    print(f"输出: {OUT_DIR}")

    # ─── Output SRA download commands ───
    sra_cmds = []
    for ds in sra_datasets:
        sra_id = ds['sra'][0] if ds.get('sra') else ''
        if sra_id:
            title = ds.get('title','')[:60]
            stresses = '+'.join(ds.get('stresses', []))
            sra_cmds.append(f"# [{stresses}] {title}")
            sra_cmds.append(f"prefetch {sra_id} && fastq-dump --split-files --gzip {sra_id}")
            sra_cmds.append("")

    sra_script = os.path.join(OUT_DIR, "download_sra.sh")
    with open(sra_script, 'w') as f:
        f.write("#!/bin/bash\n# SRA数据下载脚本 (共{}个数据集)\n".format(len(sra_datasets)))
        f.write("# 需要: SRA Toolkit (prefetch + fastq-dump)\n")
        f.write("# 注意: SRA数据为原始FASTQ，单个数据集可能10-50GB\n\n")
        f.write('\n'.join(sra_cmds))

    # ─── Download SRA metadata (not raw FASTQ - too large) ───
    print(f"\n━━━ SRA元数据下载 ━━━")
    sra_success = 0
    for i, ds in enumerate(sra_datasets):
        sra_id = ds['sra'][0] if ds.get('sra') else ''
        if not sra_id: continue

        out_subdir = os.path.join(OUT_DIR, sra_id)
        if os.path.exists(out_subdir) and os.listdir(out_subdir):
            sra_success += 1
            continue

        os.makedirs(out_subdir, exist_ok=True)

        # Download run metadata via ENA API (tiny XML, always works)
        ena_url = f"https://www.ebi.ac.uk/ena/browser/api/xml/{sra_id}"
        out_xml = os.path.join(out_subdir, f"{sra_id}_meta.xml")

        if i % 20 == 0:
            print(f"  [{i+1}/{len(sra_datasets)}] {sra_id} ...")

        if try_download(ena_url, out_xml):
            sra_success += 1

        time.sleep(0.3)  # Be nice to EBI

    print(f"SRA元数据下载: {sra_success}/{len(sra_datasets)}")
    print(f"\n━━━ 全部完成 ━━━")
    print(f"GEO表达数据: {success}个")
    print(f"SRA元数据: {sra_success}个")

if __name__ == '__main__':
    main()
