import requests

url = "https://finviz.com/quote.ashx?t=FB"
page = requests.get(url)
contents = page.content

print contents
