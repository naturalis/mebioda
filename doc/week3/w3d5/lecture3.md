The Standard Cross-Cultural Sample (SCCS) describes a subset of 186 societies that also 
appear in the EA. These 186 societies were chosen by George P. Murdock, the creator of 
the dataset, with the goal of maximizing independence among societies. The societies are 
globally distributed. Since the first installment of the SCCS was published by Murdock 
and Morrow in 1970, many others have followed, each typically focused on a particular 
theme (e.g., childhood, kinship systems). Over 2000 cultural variables have now been 
coded for the 186 societies in the SCCS. Of these, the following variables record staple
crop use:

- `SCCS1123` - [Major Agricultural Staple][SCCS1123]
- `SCCS1125` - [Second Agricultural Staple][SCCS1125]
- `SCCS1349` - [Primary Crop Name][SCCS1349]
- `SCCS1350` - [Other Primary Crop Name][SCCS1350]
- `SCCS1445` - [Secondary Crop Name][SCCS1445]
- `SCCS1446` - [Other Secondary Crop Name][SCCS1446]

We can extract these from the sample using the [dplace](https://github.com/D-PLACE/dplace-data)
command line interface:

```bash
dplace extract \
	--dataset SCCS \
	--variable SCCS1123,SCCS1125,SCCS1349,SCCS1350,SCCS1445,SCCS1446 \
	lecture3/crop_usage.csv
```

This gives us a large [table](lecture3/crop_usage.csv) with the following columns:

- **ID** - primary key, useful for dictionaries, e.g. `SCCS1`
- XD_ID
- Glottocode
- **Name** - short name for display, e.g. `Nama` - note: UTF-8
- OriginalName
- FocalYear
- Latitude
- Longitude
- **Variable** - one of the above variables, e.g. `SCCS1123`
- **Value** - one of the values, e.g. from [SCCS1123](lecture3/SCCS1123.csv)

We want to transform this to a table where each row is a society (ID) and the columns
are the different crops. Here's a [script](lecture3/transform.py) to do this in python,
which results in the following [table](lecture3/checklist.csv)

- Remove non-crop columns, e.g. `Agriculture not practiced or confined to non-food crops`
- Remove non-staple crops, e.g. `Cardamum`
- Merge alternate spellings, e.g. `Wet Rice` and `Wet rice`

Which results in this [cleaned table](lecture3/checklist-cleaned.csv)


----
[SCCS1123]: lecture3/SCCS-variables.md#SCCS1123
[SCCS1125]: lecture3/SCCS-variables.md#SCCS1125
[SCCS1349]: lecture3/SCCS-variables.md#SCCS1349
[SCCS1350]: lecture3/SCCS-variables.md#SCCS1350
[SCCS1445]: lecture3/SCCS-variables.md#SCCS1445
[SCCS1446]: lecture3/SCCS-variables.md#SCCS1446