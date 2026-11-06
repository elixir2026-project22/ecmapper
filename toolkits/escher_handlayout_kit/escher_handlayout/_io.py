"""Tiny file helpers that always close their files."""
import json

def read_text(path, mode='r'):
    with open(path, mode, **({} if 'b' in mode else {'encoding': 'utf-8'})) as f: return f.read()

def write_text(path, text):
    with open(path, 'w', encoding='utf-8') as f: f.write(text)

def read_json(path):
    with open(path, encoding='utf-8') as f: return json.load(f)

def write_json(path, obj, indent=1):
    with open(path, 'w', encoding='utf-8') as f: json.dump(obj, f, indent=indent)
