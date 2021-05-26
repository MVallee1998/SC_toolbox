import sqlite3
import SimplicialComplex as sc
import json
import datetime

m = 14
n = 10

def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data

# setting up the connection with database
conn = sqlite3.connect("Real_Toric_Spaces.db")
print("Successfully connected to the database")
cursor = conn.cursor()
facets_list = read_file('./partial_results/PLS_%d_%d' % (m, n))
# storing values in a variable
values = []
query = "INSERT OR IGNORE INTO PLSpheres (n, m, Pic, max_faces, min_non_faces, f_vector, h_vector, g_vector, " \
        "who_found_me, Date) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?); "

for k in range(len(facets_list)):
    if k%10==0:
        print((k / len(facets_list)) * 100)
    facets = facets_list[k]
    facets_converted = json.loads(facets)
    K = sc.PureSimplicialComplex(facets_converted)
    K.create_f_vector()
    K.create_g_vector()
    K.create_h_vector()
    K.compute_MNF_set()
    K.MNF_bin_to_MNF()
    min_non_faces = []
    author = 'Mathieu Vallée'
    date = datetime.datetime.now()
    values.append((K.n, K.m, K.Pic, str(K.facets), str(K.MNF_set), str(K.f_vector), str(K.h_vector), str(K.g_vector), author,
                   date.strftime("%x,%X")))
# executing the query with values

cursor.executemany(query, values)
# to make final output we have to run the 'commit()' method of the database object

conn.commit()
print(cursor.rowcount, "records inserted")
# closing the connection
cursor.close()
conn.close()