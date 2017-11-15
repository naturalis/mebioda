import unicodecsv as csv

society = {}
crops = {}
mapping = {}


# reads an individual mapping file, has 'ID' and 'Description' columns
# allows lookups of the crop name through the id
def read_mapping(variable_id):
    path = variable_id + '.csv'
    with open(path) as csvfile:
        mapping[variable_id] = {}
        mapping_reader = csv.DictReader(csvfile, delimiter=',')
        for crop_map in mapping_reader:
            crop_name = crop_map['Description']
            mapping[variable_id][crop_map['ID']] = crop_name
            if not crop_name in crops:
                crops[crop_name] = {}
            crops[crop_name][variable_id] = crop_map['ID']


# open main spreadsheet
with open('crop_usage.csv') as csvfile:

    # read rows as dictionaries
    reader = csv.DictReader(csvfile, delimiter=',')
    for row in reader:

        # society ID, usage variable, value
        soc = row['ID']
        usage = row['Variable']
        value = row['Value']

        # first time we see this society
        if not soc in society:
            society[soc] = {}

        # first time we see this usage variable
        if not usage in mapping:
            read_mapping(usage)
        crop = mapping[usage][value]
        society[soc][crop] = 1

    # print output
    crop_names = sorted(crops.keys())
    header = ['"ID']
    header.extend(crop_names)
    print '","'.join(header) + '"'
    for soc in sorted(society.keys()):
        row = [soc]
        for crop in crop_names:
            if crop in society[soc]:
                row.append("1")
            else:
                row.append("0")
        print ','.join(row)

