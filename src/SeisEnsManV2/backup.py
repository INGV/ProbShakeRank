import os
from datetime import datetime, timedelta
import json
import shutil

# BACKUP of 'event_stat.json'
file_dir = os.path.join(os.getcwd(), 'input')
eventstatjson = 'event_stat.json'
eventstatjson_fullpath = os.path.join(file_dir, eventstatjson)
print("eventstatjson_fullpath: ", eventstatjson_fullpath)

timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
if os.path.exists(eventstatjson_fullpath):
    backup_folder = os.path.join(os.path.join(os.getcwd()), 'BACKUP/EVENT_STAT')
    if not os.path.exists(backup_folder):
        os.makedirs(backup_folder)

    backup_filename = f"event_stat_{timestamp}.json"
    backup_full_path = os.path.join(backup_folder, backup_filename)

    shutil.copy(eventstatjson_fullpath, backup_full_path)

# BACKUP OF 'output' FOLDER
start_folder = os.path.join(os.getcwd(), 'output')
if not os.path.exists(start_folder):
    os.makedirs(start_folder)
timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
backup_folder = os.path.join(os.path.join(os.getcwd()), f"BACKUP/LIST_SCENARIOS_{timestamp}")
if not os.path.exists(backup_folder):
    os.makedirs(backup_folder)
for item in os.listdir(start_folder):
    item_path = os.path.join(start_folder, item)
    if os.path.isfile(item_path) and item != ".DS_Store":
        shutil.move(item_path, backup_folder)
