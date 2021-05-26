import sqlite3
import SimplicialComplex
import json
import datetime

def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data
# setting up the connection with database
conn = sqlite3.connect("BAK/Real_Toric_Spaces.db")
print("Successfully connected to the database")
cursor = conn.cursor()
facets_list = read_file('.Result_Hyuntae/PLS_10_6')
# storing values in a variable
values = []
k = 0
query2 = "INSERT OR IGNORE INTO RealToricSpaces (FK_PLS_Key, Z2_Char_function) VALUES (?, ?); "

for row in cursor.execute('SELECT PLS_Key,max_faces FROM PLSpheres'):
    print(k / 6713)
    k += 1
    K_Key = row[0]
    K = json.loads(row[1])
    char_funct_over_K = SimplicialComplex.Garrison_Scott(K)
    if len(char_funct_over_K) != 0:
        for char_funct in char_funct_over_K:
            values.append((K_Key, str(char_funct)))
cursor.executemany(query2, values)
conn.commit()
print(cursor.rowcount, "records inserted")