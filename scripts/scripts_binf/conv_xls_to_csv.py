import pandas as pd
import sys 
import os 
if len(sys.argv) != 3:
    print("Usage: conv_xls_to_csv.py <input.xls> <output.csv>" )
    sys.exit(1)  

def conv_xls_to_csv():
   """
   Converts .xls-file to .csv-file
   """
   #read an excel file and convert into a dataframe object
   df = pd.DataFrame(pd.read_excel("{0}".format(sys.argv[1])))
  
   # convert to csv
   return df.to_csv(sys.argv[2])

if __name__ == '__main__':
   conv_xls_to_csv()


