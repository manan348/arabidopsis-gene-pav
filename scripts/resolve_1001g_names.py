"""Map 1001 Genomes numeric IDs to accession names via AraPheno REST."""
import json, sys, urllib.request
from concurrent.futures import ThreadPoolExecutor

def fetch(i):
    try:
        with urllib.request.urlopen(f"https://arapheno.1001genomes.org/rest/accession/{i}.json", timeout=30) as r:
            d = json.load(r)
        return f"{i}\t{d['name']}\t{d['country']}"
    except Exception:
        return f"{i}\tNA\tNA"

ids = [l.strip() for l in open(sys.argv[1]) if l.strip()]
with ThreadPoolExecutor(16) as ex, open(sys.argv[2], "w", encoding="utf-8") as out:
    for line in ex.map(fetch, ids):
        out.write(line + "\n")
