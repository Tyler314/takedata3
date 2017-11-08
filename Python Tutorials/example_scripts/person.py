
class Person:
    def __init__(self, name, job=None, pay=0):
        self.name = name
        self.job = job
        self.pay = pay

    def lastName(self):
        return self.name.split()[-1]

    def giveRaise(self, percent):
        self.pay = round(self.pay * (1 + percent), 2)

    # Operator Overloading
    def __gt__(self, other):
        return self.pay > other.pay

    def __lt__(self, other):
        return self.pay < other.pay

    def __eq__(self, other):
        return self.pay == other.pay

    def __add__(self, other):
        return Person(name = self.name + ' and ' + other.name, pay = self.pay + other.pay)