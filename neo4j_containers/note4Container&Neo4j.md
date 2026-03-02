# Notes of using docker and neo4j

### Multiple container
- docker-compose.yml file: can write multiple graph container
    1. Port Mapping (Host:Container)
        + Containers live in their own network bubble.
        + For the first container: 7474:7474 means you go to localhost:7474.
        + For the second container: 7575:7474 means you go to localhost:7575.
        + Important: Inside the container, Neo4j still thinks it's on 7474, but Docker acts as a translator.
    2. Unique Volume Names
        + If both containers used the volume neo4j:/data, they would try to write to the same database files, causing a crash or data corruption. I renamed them to neo4j_data_organic and neo4j_data_kg.
    3. Container Names
        + The container_name property ensures that when you run docker ps, you can easily tell which graph is which.
    4. Execuation:
       + in terminal run `docker-compose up -d` can start all the defined containers
       + open localhost:7474 (or other port) can open neo4j browser

### Upload KG to neo4j-database
#### Cypher file
+ run `cat file_name.cypher | docker exec -i container_name cypher-shell -u neo4j
 -p test1234` in terminal

#### pkl file
- connect to localhost by a driver
- load KG
- write the triples to neo4j-database