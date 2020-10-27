class Node:

    def __init__(self, data=None):

        self.left = None
        self.right = None
        self.data = data

    def isEmpty(self):
        return self.data == None

    def size(self):
        if self.isEmpty():
            return 0
        left_size = 0
        right_size = 0
        if self.left:
            left_size = self.left.size()
        if self.right:
            right_size = self.right.size()
        return 1 + left_size + right_size

    def insert(self, data):
        # Compare the new value with the parent node
        if self.data:
            if data < self.data[0]:
                if self.left is None:
                    self.left = Node(data)
                else:
                    self.left.insert(data)
            elif data > self.data[0]:
                if self.right is None:
                    self.right = Node(data)
                else:
                    self.right.insert(data)
            elif data == self.data[0]:
                self.data[1] += data[1]
        else:
            self.data = data

    # insert a list of values
    def insertList(self, l):
        for val in l:
            self.insert(val)

    # findval method to compare the value with nodes
    def findval(self, lkpval):
        if self.isEmpty():
            return False
        if lkpval < self.data[0]:
            if self.left is None:
                return False
            return self.left.findval(lkpval)
        elif lkpval > self.data[0]:
            if self.right is None:
                return False
            return self.right.findval(lkpval)
        else:
            return True

    # Print the tree
    def PrintTree(self):
        if self.left:
            self.left.PrintTree()
        print(self.data),
        if self.right:
            self.right.PrintTree()

    # Transform the tree into a list
    def TreeToList(self, l):
        if self.left:
            self.left.TreeToList(l)
        if self.data:
            l.append(self.data),
        if self.right:
            self.right.TreeToList(l)
