from person import Person

class GradStudent(Person):
    # Static Variable
    numInstances = 0
    def __init__(self, name, job='TA', pay=0, areaOfStudy = 'Plasma'):
        Person.__init__(self, name, job, pay)
        # Inheritance
        self.areaOfStudy = areaOfStudy
        GradStudent.makeInstance()

    def timeLeft(self):
        if self.areaOfStudy in ['Plasma', 'Nuclear']:
            base = 9
        elif self.areaOfStudy in ['High Energy', 'Atomic']:
            base = 7
        else:
            base = 5
        reduce = 2 if self.job == 'TA' else 4
        return base - reduce

    # Polymorphism
    def giveRaise(self, percent):
        return Person.giveRaise(self, percent * 0.9)

    @staticmethod
    def makeInstance():
        GradStudent.numInstances += 1
        print('Created {} instances of GradStudent'.format(GradStudent.numInstances))
        