import boto3
import os

def download_all_refs(bucket="",organism="Staphylococcus_aureus", outdir="bb_cluster_references"):
    """
    Download all ref.fa.gz files from s3://bucket/dbs/<organism>/<cluster>/ref.fa.gz
    and save them locally as <organism>_<cluster>_ref.fa.gz
    """
    s3 = boto3.client("s3")
    prefix = f"waphl-nextflow-batch/groups/general/bigbacter/db/{organism}/clusters"
    os.makedirs(outdir, exist_ok=True)

    paginator = s3.get_paginator("list_objects_v2")
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        if "Contents" not in page:
            continue
        for obj in page["Contents"]:
            key = obj["Key"]
            #print(obj)
            if key.endswith("ref.fa.gz"):
                # Expect keys like: dbs/<organism>/<cluster>/ref.fa.gz
                parts = key.split("/")
                if len(parts) >= 4:
                    
                    cluster = [parts[i + 1] for i in range(len(parts) - 1) if parts[i] == "clusters"][0]
                    print(cluster,"cluster", key)
                    local_name = f"{organism}_{cluster}_ref.fa.gz"
                    local_path = os.path.join(outdir, local_name)
                    print(f"Downloading {key} → {local_path}")
                    s3.download_file(bucket, key, local_path)

# Example usage:
# download_all_refs("my-bucket")
if __name__ == "__main__":
    download_all_refs("waphl-nextflow-batch-bucket",sys.argv[1])