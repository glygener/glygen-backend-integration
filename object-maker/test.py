import json
import requests

url = "https://api.glygen.org/misc/verlist"
res = requests.get(url)
res_obj = json.loads(res.content)
   
print (json.dumps(res_obj, indent=4))

