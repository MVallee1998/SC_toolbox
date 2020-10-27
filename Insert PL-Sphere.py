import sqlite3
import SimplicialComplex
import json
import datetime


# setting up the connection with database
conn = sqlite3.connect("Real_Toric_Spaces.db")
print("Successfully connected to the database")
cursor = conn.cursor()
facets_list = SimplicialComplex.read_file('./PLS_10_6')
# storing values in a variable
values = []
query = "INSERT OR IGNORE INTO PLSpheres (n, m, Pic, max_faces, min_non_faces, f_vector, h_vector, g_vector, " \
        "who_found_me, Date) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?); "

for k in range(len(facets_list)):
    print((k / len(facets_list)) * 100)
    facets = facets_list[k]
    facets_converted = json.loads(facets)
    n, m, pic = SimplicialComplex.give_dim(facets_converted)
    faces_set = SimplicialComplex.create_faces_set(facets_converted)
    f_vector = SimplicialComplex.create_f_vector(faces_set)
    h_vector = SimplicialComplex.create_h_vector(f_vector)
    g_vector = SimplicialComplex.create_g_vector(h_vector)
    MNF_set = SimplicialComplex.faces_set_to_MNF_set(faces_set)
    min_non_faces = []
    MNF_set.TreeToList(min_non_faces)
    author = 'Hyuntae Jang'
    date = datetime.datetime.now()
    values.append((n, m, pic, facets, str(min_non_faces), str(f_vector), str(h_vector), str(g_vector), author,
                   date.strftime("%x,%X")))
# executing the query with values

cursor.executemany(query, values)
# to make final output we have to run the 'commit()' method of the database object

conn.commit()
print(cursor.rowcount, "records inserted")
# closing the connection
cursor.close()
conn.close()