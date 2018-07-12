import urllib
import simplejson as json # sudo pip install simplejson
url = "http://www.boldsystems.org/index.php/API_Tax/TaxonSearch?taxName=Artiodactyla"
response = urllib.urlopen(url)
data = json.loads(response.read())

if data['top_matched_names']:
	for name in data['top_matched_names']:
		if name['representitive_image']:
			print name['representitive_image']['image']