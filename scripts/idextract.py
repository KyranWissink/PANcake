import json

def extract_node_ids(json_file):
    node_ids = []
    with open(json_file, 'r') as file:
        for line in file:
            alignment = json.loads(line)
            if 'path' in alignment and 'mapping' in alignment['path']:
                for mapping in alignment['path']['mapping']:
                    if 'position' in mapping and 'node_id' in mapping['position']:
                        node_id = mapping['position']['node_id']
                        node_ids.append(node_id)
    return node_ids

# Replace 'giraffe_4287_effectors.json' with the path to your JSON file
node_ids = extract_node_ids('annotated.json')

# Print the extracted node IDs
for nid in set(node_ids):  # Using set to remove duplicate IDs
    print(nid)

