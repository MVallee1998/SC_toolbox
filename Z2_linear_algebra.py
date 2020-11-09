class Z2Array:
    def __init__(self, m, values):
        self.n = len(values)
        self.m = m
        self.values = values

    def is_invertible(self):
        if self.m != self.n:
            raise Exception
        else:
            list_2_pow = [1]
            for k in range(1, self.n):
                list_2_pow.append(list_2_pow[-1] * 2)
            for i in range(0, self.n):
                j_0 = -1
                for j in range(i, self.n):
                    if (self.values[j] | list_2_pow[i]) == self.values[j]:
                        j_0 = j
                        break
                if j_0 == -1:
                    return False
                else:
                    (self.values[i], self.values[j_0]) = (self.values[j_0], self.values[i])
                for j in range(i + 1, self.n):
                    if (self.values[j] | list_2_pow[i]) == self.values[j]:
                        self.values[j] = self.values[j] ^ self.values[i]
            for a in self.values:
                if a == 0:
                    return False
            return True

    def Z2_rank(self):
        list_2_pow = [1]
        for k in range(1, self.m):
            list_2_pow.append(list_2_pow[-1] * 2)

        i = 0
        j = 0
        while i < self.n and j < self.m:
            j_0 = -1
            for j_ind in range(i, self.m):
                if (self.values[i] | list_2_pow[j_ind]) == self.values[i]:
                    j_0 = j_ind
                    break
            if j_0 == -1:
                i += 1
                continue
            else:
                for l in range(self.n):
                    if (self.values[l] | list_2_pow[j_0] == self.values[l]) ^ (self.values[l]| list_2_pow[i] == self.values[l]):
                        self.values[l] ^= list_2_pow[i]
                        self.values[l] ^= list_2_pow[j_0]

                for j_ind in range(i + 1, self.m):
                    if (self.values[i] | list_2_pow[j_ind]) == self.values[i]:
                        for l in range(self.n):
                            if (self.values[l] | list_2_pow[i] == self.values[l]):
                                self.values[l] ^= list_2_pow[i]
                i += 1
                j += 1
        return j
