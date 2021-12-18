class BTS:
    def __init__(self, ind, val):
        self.left = None
        self.right = None
        self.index = ind
        self.value = val

    def insert(self, ind, val):
        if self.value:
            if val < self.value:
                if self.left is None:
                    self.left = BTS(ind, val)
                else:
                    self.left.insert(ind, val)
            elif val > self.value:
                if self.right is None:
                    self.right = BTS(ind, val)
                else:
                    self.right.insert(ind, val)
            else:
                self.index = ind
                self.value = val

    def search(self, val):
        if self.value:
            if val < self.value:
                if self.left:
                    return self.left.search(val)
                else:
                    return -1
            elif val > self.value:
                if self.right:
                    return self.right.search(val)
                else:
                    return -1
            else:
                return self.index
