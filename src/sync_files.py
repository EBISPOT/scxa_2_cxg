import os

# Define the remote host and base directories
remote_host = "FILL_ME_IN"
remote_base_dir = "FILL_ME_IN"
download_dir = os.path.join("downloads")


exp_list = ["E-MTAB-8698",
            "E-MTAB-9444",
            "E-MTAB-10519",
            "E-MTAB-10628",
            "E-MTAB-7194",
            "E-MTAB-7195",
            "E-GEOD-141273",
            "E-GEOD-100058",
            "E-GEOD-134722",
            "E-GEOD-136162",
            "E-GEOD-146040",
            "E-GEOD-172231",
            "E-GEOD-152495",
            "E-GEOD-141807",
            "E-GEOD-157775",
            "E-GEOD-103771",
            "E-GEOD-126139",
            "E-GEOD-157202",
            "E-GEOD-125948",
            "E-GEOD-125948",
            "E-GEOD-147601",
            "E-CURD-21",
            "E-CURD-134",
            "E-CURD-91",
            "E-CURD-90",
            "E-CURD-92",
            "E-CURD-56",
            "E-CURD-124",
            "E-CURD-87"]

# Iterate over the list of experiment IDs
for exp_id in exp_list:
    # Construct the source and destination paths
    # remote_path = f"{remote_base_dir}"
    # local_path = os.path.join(download_dir, exp_id)
    # filename = f"{local_path}/{exp_id}_modified.project.h5ad"
    # filename = f"{local_path}/{exp_id}.owl"

    # Construct the rsync command
    # rsync_command = f"rsync -avz --progress {filename} {remote_host}:{remote_base_dir}"
    # cp_command = f"cp {filename} graphs"
    # Execute the rsync command
    print(f"https://github.com/EBISPOT/scxa_2_cxg/raw/main/graphs/{exp_id}.owl")
    # os.system(rsync_command)
    # os.system(cp_command)
    # print(filename)
    # print(os.path.exists(filename))

print("Synchronization complete.")
