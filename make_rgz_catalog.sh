#!/bin/bash

RGZ_PATH='/data/tabernacle/larry/RGZdata'

# Start MongoDB
echo "Starting MongoDB"
numactl --interleave=all mongod --fork --logpath $RGZ_PATH"/mongodb/log/mongodb.log" --dbpath $RGZ_PATH"/data/db"

# Restore the raw files
echo "Restoring MongoDB files"

arr=($(find $RGZ_PATH"/mongodb/exports/sanitized_radio"* -type d))
BACKUP_PATH=${arr[${#arr[@]}-1]}
echo ${#arr[@]}" backups of catalog found"
echo "Using "$BACKUP_PATH

#mongorestore --db radio --drop --collection radio_users $BACKUP_PATH'/radio_users.bson'

mongoimport --db radio --drop --collection radio_subjects $BACKUP_PATH'/radio_subjects.json'
mongoimport --db radio --drop --collection radio_classifications $BACKUP_PATH'/radio_classifications.json'
mongoimport --db radio --drop --collection radio_groups $BACKUP_PATH'/radio_groups.json'

# Restore the latest version of the catalog
mongorestore --db radio --drop --collection catalog $RGZ_PATH'/rgz_mongo/radio/catalog.bson'

# Activate necessary Python environments

echo "Activating Python environments"
module load python27
source /home/phys/willett/veastropy/bin/activate

# Run consensus algorithm

echo "Running consensus algorithm on new subjects"
python2.7 $RGZ_PATH"/rgz-analysis/python/consensus.py"

# Update catalog

echo "Updating catalog"
python2.7 $RGZ_PATH"/rgz-analysis/RGZcatalogCode/RGZcatalog.py" --consensus $RGZ_PATH"/rgz-analysis/csv/consensus_rgz_first.csv"

# Create static version (full flat catalog, removing peak/component data)

echo "Outputting static catalog"
python2.7 $RGZ_PATH"/rgz-analysis/python/static_catalog.py"

# Close Python environment

echo "Closing Python environments"
deactivate
module unload python27

# Export the Mongo files

echo "Dumping Mongo files"
mongodump --db radio --out rgz_mongo

# Kill the MongoDB process

echo "Killing MongoDB"
kill $(ps aux | grep '[m]ongod --fork' | awk '{print $2}')
