#!/bin/bash

## Update MongoDB with the latest data from Radio Galaxy Zoo. 
## Must supply the name of the directory where the BSON files are stored as input

# Example:
# >> ./update_rgz_mongo.sh radio_2014-04-03/

echo 'Starting mongod'

mongod --fork --logpath log/mongodb.log

echo 'Restoring files'
mongorestore --db ouroboros --drop --collection radio_users $1'radio_users.bson'
mongorestore --db ouroboros --drop --collection radio_subjects $1'radio_subjects.bson'
mongorestore --db ouroboros --drop --collection radio_classifications $1'radio_classifications.bson'

echo 'Killing the mongod process'

kill $(ps aux | grep '[m]ongod --fork' | awk '{print $2}')
