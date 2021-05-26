import mysql.connector as mysql
from datetime import date, datetime, timedelta


def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data


db = mysql.connect(
    host="localhost",
    user="root",
    passwd="",
    database="Real Toric Manifolds"
)

cursor = db.cursor()

query = "INSERT INTO pure_simplicial_complex (n, m, Pic, max_faces, is_PL_Sphere, who_found_me)" \
        "VALUES (%s, %s, %s, %s, %s, %s)"

# storing values in a variable
max_faces_list = read_file('./PLS_9_5.txt')
values = []
for max_face in max_faces_list:
    values.append((5, 9, 4, max_face, 1, 'Hyuntae'))

# executing the query with values
cursor.executemany(query, values)

# to make final output we have to run the 'commit()' method of the database object
db.commit()

print(cursor.rowcount, "records inserted")
