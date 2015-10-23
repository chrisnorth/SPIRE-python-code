from numpy import zeros
from astropy.table import Table,Column
from matplotlib import is_string_like

tabIn = Table.read('stats/Results-floats.csv', format='ascii.csv')
tabOut=Table()
for col in tabIn.colnames:
    if is_string_like(tabIn[col][0]):
        tabOut.add_column(Column(zeros(len(tabIn)),name=col))
        for i in range(len(tabIn)):
            tabOut[col][i]=float(tabIn[col][i].replace(',',''))
    else:
        tabOut.add_column(tabIn[col])
tabOut.write('stats/Results-floats2.csv', format='ascii.csv')
