import sqlite3
import SimplicialComplex as sc
import json
import numpy as np
import datetime
m = 11
n = 7

list_2_pow = np.ones(16,dtype=int)
for k in range(15):
    list_2_pow[k+1] = list_2_pow[k]*2


def read_file(filename):
    with open(filename, 'rb') as f:
        data = f.readlines()
    data = [x.strip() for x in data]
    return data

# setting up the connection with database
conn = sqlite3.connect("Real_Toric_Spaces.db")
print("Successfully connected to the database")
cursor = conn.cursor()
facets_list = read_file('./final_results_BAK/PLS_%d_%d' % (m, n))

# storing values in a variable
values = []
counter = 0
query = "INSERT OR IGNORE INTO RealToricSpaces (FK_PLS_Key, Z2_Char_function) VALUES (?, ?); "
for row in cursor.execute('SELECT PLS_Key,max_faces FROM PLSpheres'):
    counter += 1
    if counter%100==0:
        print(counter/8794)
    K_Key = row[0]
    K_MF = json.loads(row[1])
    K = sc.PureSimplicialComplex(K_MF)
    if True:
        p = K.Pic
        n = K.n
        m = K.m
        list_IDCM_bin = sc.IDCM_Garrison_Scott(K)
        list_IDCM= []
        for IDCM_bin in list_IDCM_bin:
            IDCM = np.zeros((m,p))
            for i in range(m):
                for j in range(p):
                    if list_2_pow[j] | IDCM_bin[i] == IDCM_bin[i]:
                        IDCM[i,j] = 1
            list_IDCM.append(IDCM.copy())
        list_CM= []
        for IDCM in list_IDCM:
            CM = np.zeros((n,m))
            CM[:,:n] = np.eye(n)
            CM[:,n:m] = IDCM[:n,:]
            list_CM.append(CM.copy())
        list_CM_bin = []
        for CM in list_CM:
            CM_bin = []
            for k in range(m):
                CM_bin.append(np.sum(list_2_pow[np.flatnonzero(CM[:,k])]))
            list_CM_bin.append(CM_bin.copy())
        if len(list_CM_bin) != 0:
            for char_funct in list_CM_bin:
                values.append((K_Key, str(char_funct)))
print(counter)
cursor.executemany(query, values)
conn.commit()
print(cursor.rowcount, "records inserted")
conn.close()