"""Portable read-only streaming queries; all result fields remain literal strings."""
import argparse
import csv
import gzip
import hashlib
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent

def load():
    manifest=json.loads((ROOT/'manifest.json').read_text(encoding='utf-8'))
    for table in manifest['tables']:
        expected_start=1
        for physical in physical_files(table):
            path=(ROOT/physical['path']).resolve()
            if not path.is_relative_to(ROOT) or pathlib.PurePosixPath(physical['path']).is_absolute():
                raise ValueError('Manifest path escapes repository')
            if 'parts' in table:
                if physical['source_row_start_1based']!=expected_start or physical['row_count']<1 or physical['source_row_end_1based']!=expected_start+physical['row_count']-1:
                    raise ValueError('Part row ranges are not contiguous')
                expected_start=physical['source_row_end_1based']+1
        if 'parts' in table and expected_start!=table['published_rows']+1:
            raise ValueError('Logical row range mismatch')
    return manifest

def physical_files(table):
    return table['parts'] if 'parts' in table else [table]

def rows(table):
    for physical in physical_files(table):
        start=physical.get('source_row_start_1based',1)
        with gzip.open(ROOT/physical['path'],'rt',encoding='utf-8-sig',newline='') as f:
            yield from enumerate(csv.DictReader(f,delimiter=table['delimiter']),start)

def match(table,row,mode,value):
    roles={'feature':{'feature_identifier','source_identifier'},'gene':{'symbol','source_gene_label'},'pathway':{'pathway_identifier','pathway_name'}}[mode]
    for col,field in table['fields'].items():
        role=field.get('role')
        if mode=='pathway':
            eligible=role in roles or col in {'pathway','term','term_name'}
        else: eligible=role in roles
        if not eligible: continue
        literal=row.get(col,'')
        values=literal.split(field['delimiter']) if field.get('delimiter') else [literal]
        if mode=='feature':
            if literal==value: return True
        elif any(x.strip().casefold()==value.casefold() for x in values): return True
    return False

def validate(manifest):
    checked=0; row_total=0
    for table in manifest['tables']:
        total=0; reconstructed=hashlib.sha256()
        for index,physical in enumerate(physical_files(table)):
            path=ROOT/physical['path']; h=hashlib.sha256()
            with path.open('rb') as f:
                for b in iter(lambda:f.read(1048576),b''): h.update(b)
            if h.hexdigest()!=physical['published_sha256']: raise ValueError('File SHA mismatch: '+table['table_id'])
            count=0
            with gzip.open(path,'rt',encoding='utf-8-sig',newline='') as f:
                reader=csv.reader(f,delimiter=table['delimiter'])
                if next(reader)!=table['columns']: raise ValueError('Header mismatch')
                for count,row in enumerate(reader,1):
                    if len(row)!=len(table['columns']): raise ValueError('Row width mismatch')
            if count!=physical.get('row_count',table['published_rows']): raise ValueError('Row count mismatch')
            total+=count
            if 'parts' in table:
                with gzip.open(path,'rb') as f:
                    header=f.readline()
                    if index==0: reconstructed.update(header)
                    for b in iter(lambda:f.read(1048576),b''): reconstructed.update(b)
        if total!=table['published_rows']: raise ValueError('Logical row count mismatch')
        if 'parts' in table and reconstructed.hexdigest()!=table['logical_reassembled_decompressed_sha256']: raise ValueError('Logical reconstruction hash mismatch')
        checked+=1; row_total+=total
    return {'status':'PASS','checked_tables':checked,'rows':row_total,'scientific_complete':False}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['status','tables','feature','gene','pathway','validate'])
    parser.add_argument('value',nargs='?')
    parser.add_argument('--package'); parser.add_argument('--table')
    parser.add_argument('--limit',type=int,default=100); parser.add_argument('--offset',type=int,default=0)
    parser.add_argument('--all',action='store_true')
    args=parser.parse_args()
    if args.limit<1 or args.offset<0: parser.error('limit must be positive and offset nonnegative')
    if args.mode in {'feature','gene','pathway'} and not args.value: parser.error('query value required')
    if args.mode not in {'feature','gene','pathway'} and args.value: parser.error('unexpected query value')
    manifest=load()
    if args.mode=='validate': print(json.dumps(validate(manifest))); return
    if args.mode=='status':
        print(json.dumps({k:v for k,v in manifest.items() if k!='tables'})); return
    tables=[t for t in manifest['tables'] if (not args.package or t['package']==args.package) and (not args.table or t['table_id']==args.table)]
    if args.mode=='tables':
        values=[{k:v for k,v in t.items() if k not in {'fields'}} for t in tables]
        print(json.dumps(values[args.offset:] if args.all else values[args.offset:args.offset+args.limit])); return
    aliases={'PD-1':'PDCD1','PD-L1':'CD274'}
    resolved=aliases.get(args.value.upper(),args.value) if args.mode=='gene' else args.value
    emitted=0; skipped=0
    for table in tables:
        for number,row in rows(table):
            if not match(table,row,args.mode,resolved): continue
            if skipped<args.offset: skipped+=1; continue
            print(json.dumps({'table_id':table['table_id'],'package':table['package'],'source_row_1based':number,'query':args.value,'resolved_query':resolved,'alias_mapping':resolved!=args.value,'status':table['status'],'record_role':table['record_role'],'fields':row},ensure_ascii=False))
            emitted+=1
            if not args.all and emitted>=args.limit: return

if __name__=='__main__':
    try: main()
    except (OSError,ValueError,KeyError) as error:
        print(str(error),file=sys.stderr); sys.exit(1)
