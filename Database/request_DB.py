import sqlite3
import SimplicialComplex
import json

def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data
# setting up the connection with database
conn = sqlite3.connect("BAK/Real_Toric_Spaces.db")
print("Successfully connected to the database")
cursor = conn.cursor()
# facets_list = read_file('.Result_Hyuntae/CSPLS_6_10')
# storing values in a variable
values = []
k = 0
# query = "INSERT OR IGNORE INTO RealToricSpaces (FK_PLS_Key, Z2_Char_function) VALUES (?, ?); "

for row in cursor.execute('SELECT PLS_Key, m, min_non_faces FROM PLSpheres'):
    K_Key = row[0]
    m = row[1]
    mnf = json.loads(row[2])
    if len(mnf)>m and m<8:
        print(m,mnf)

