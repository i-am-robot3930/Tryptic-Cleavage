from pylab import *
from Bio import SeqIO
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import pandas as pd
import re
import numpy as np
pd.options.mode.chained_assignment = None 
import io
from io import StringIO
import matplotlib.pyplot as plt
import seaborn as sns 
import gzip
import shutil
import glob
import scipy.stats as stats
import pylab as pl
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.metrics import classification_report
from sklearn.tree import plot_tree
import shap
import os
from Bio.PDB.ResidueDepth import ResidueDepth
from Bio.PDB.SASA import ShrakeRupley
import csv




def list_to_string(list): #convert list to string  
    return ','.join(list)

def check_Sequence_in_Protein_Sequence(row): #function to check if sequence is in full sequence
    if row['Sequence'] in row['Protein_Sequence']:
        return True
    else:
        return False


def peptide_start_position(protein, peptide): #This function finds the start position of the peptide in the protein sequence   
    for p in peptide:
        start_position= protein.find(peptide) #find start position of peptide in protein sequence  
        return start_position+1 #return start position of peptide in protein sequence

def peptide_end_position(protein, peptide): #This function finds the end position of the peptide in the protein sequence
    for p in peptide: 
        start_position= protein.find(peptide) #find start position of peptide in protein sequence
        end_position= start_position+len(peptide)-1 #find end position of peptide in protein sequence
        return end_position+1 #return end position of peptide

def find_peptide(protein, peptide): #This function finds the next AA after the peptide in the protein sequence
    try:
        position = protein.find(peptide) #find position of peptide in protein sequence
        next_position= position+len(str(peptide)) #find next position of peptide in protein sequence
        return protein[next_position]   #return next AA of peptide
        
    except:
        return ''

def find_last_AA(protein, peptide): #This function finds the last AA of the peptide in the protein sequence
    
    try:
        position = protein.find(peptide) #find position of peptide in protein sequence
        last_position_amino_acid= position + len(peptide)-1 #find last position of peptide in protein sequence

        return protein[last_position_amino_acid] #return last AA of peptide
    except:
        return ''

def next_2_AA(protein, peptide): 
    "This function finds the next 2 AA after the peptide in the protein sequence"
    try:
        position = protein.find(peptide) #find position of peptide in protein sequence
        last_position_amino_acid= position + len(peptide)-1 #find last position of peptide in protein sequence
        next_position= position+len(str(peptide)) #find next position of peptide in protein sequence
        return protein[last_position_amino_acid:next_position+1] #return next 5 AA of peptide
    except:
        return ''
    

def split_Cleavage_sites(row):  #function to split  cleavage sites
    if row == 'nan':
        return row
    else:
        return row.split(',')

def second_split(cs_split): #function to seperate 2nd cleavage site
    if len(cs_split) >2:
        return cs_split[1]
    else:
        return ''
    
def cs_split(first_half, start, Missed_Cleavage_AA):
    if Missed_Cleavage_AA == '':
        return 0 
    else:
        return start + len(first_half)-1
        
def split_second(first_half, start, Missed_Cleavage_AA):
    if Missed_Cleavage_AA == '':
        return 0 
    else:
        return start + len(first_half)
    
def fill_empty_MS_AA(row):
    if row['MS_AA'] == '':
        last_amino_acid = row['Sequence'][-1]
        return last_amino_acid
    else:
        return row['MS_AA']
    
    

def change_values(row): #function to change values
    if row['duplicate'] == True:
        row['missed'] = 0
        row['MS_AA'] = row['Sequence'][-1]
        row['position_MS'] = row['end']
        return row
    else:
        return row
        
def fill_empty_position_MS(row): #function to fill empty position_MS
    if row['position_MS'] == 0:
        row['position_MS'] = row['end']
        return row
    else:
        return row 
    


def dssp_file_df(file): # function to read in the dssp file and return a dataframe
    p = PDBParser() # create a PDB parser
    structure = p.get_structure(f"{file}", file) # get the structure from the PDB file
    model = structure[0]#  get the first model from the structure
    dssp = DSSP(model, file) # get the DSSP from the model
    df= pd.DataFrame(dssp)# create a dataframe from the DSSP
    df.columns = ["dssp index", "amino acid", "secondary structure", "relative ASA", "phi", "psi", "NH01R", "NH01E", "ONH1", "ONHE", "NH2R",
    "NO2E", "ONH2R", "ONH2E"] # naming the columns
    return df # return the dataframe

def dssp_dataframe_calculation(df):
    df.fillna('', inplace=True)# Replace all NaN values in the DataFrame with empty strings
    output_rows = [] # Create an empty list to store the output rows
    # Define the path to the proteome directory
    proteome_path = 'C:\\Users\\ahmad\\Dropbox\\PC\\Documents\\QMUL\\Research Project\\Proteomics\\Proteomics-extended\\data\\proteome' 
    
    # Iterate over each row of the input DataFrame
    for i, j in df.iterrows():
        # Generate the file path for the DSSP file for the current protein
        filename = f'AF-{j["ProteinID"]}-F1-model_v4.pdb'
        filepath = os.path.join(proteome_path, filename, filename)
        
        # Check if the DSSP file exists
        if os.path.exists(filepath):
            # If it does, read the DSSP file into a DataFrame
            df_dssp = dssp_file_df(filepath)
            # Find the row in the DSSP DataFrame that corresponds to the last amino acid in the protein sequence
            dssp_row = df_dssp[df_dssp['dssp index'] == j['LastAA_position']]
            # If a matching row was found, append the values to the output rows list
            if not dssp_row.empty:
                output_rows.append([
                    dssp_row['relative ASA'].values[0],
                    dssp_row['secondary structure'].values[0],
                    dssp_row['phi'].values[0],
                    dssp_row['psi'].values[0],
                    dssp_row['NH01R'].values[0],
                    dssp_row['NH01E'].values[0],
                    dssp_row['ONH1'].values[0],
                    dssp_row['ONHE'].values[0],
                    dssp_row['NH2R'].values[0],
                    dssp_row['NO2E'].values[0],
                    dssp_row['ONH2R'].values[0],
                    dssp_row['ONH2E'].values[0],
                ])
            else: # If no matching row was found, append 12 empty strings to the output rows list
                output_rows.append(('' for _ in range(12)))
        else:# If the DSSP file does not exist, append 12 empty strings to the output rows list
            output_rows.append(('' for _ in range(12)))
    # Write the output rows list to a CSV file named 'output.csv'
    print(len(output_rows))
    with open('output2.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['relative ASA', 'secondary_structure', 'phi', 'psi', 'NH01R', 'NH01E', 'ONH1', 'ONHE', 'NH2R', 'NO2E', 'ONH2R', 'ONH2E'])
        writer.writerows(output_rows)
    # Read the 'output.csv' file into a DataFrame
    dataframe = pd.read_csv("C:\\Users\\ahmad\\Dropbox\\PC\\Documents\\QMUL\\Research Project\\Proteomics\\Proteomics-extended\\output2.csv")
    df = df.join(dataframe)# Join the input DataFrame with the output DataFrame 
    df['relative ASA'].replace('', np.nan, inplace=True) 
    df.dropna(subset=['relative ASA'], inplace=True)
    return df


def dssp_df_calcualtion(df): # function to calculate the secondary structure and relative ASA
    relative_asa = [] # create a list to store the relative ASA
    secondary_structure = []    # create a list to store the secondary structure
    phi = []    # create a list to store the phi
    psi = []    # create a list to store the psi
    NH01R = []  # create a list to store the NH01R
    NH01E = []  # create a list to store the NH01E
    ONH1 = []   # create a list to store the ONH1
    ONHE = []   # create a list to store the ONHE
    NH2R = []   # create a list to store the NH2R
    NO2E = []   # create a list to store the NO2E
    ONH2R = []  # create a list to store the ONH2R
    ONH2E = []  # create a list to store the ONH2E



    for i,j in df.iterrows():   # loop through the dataframe
        filename= f'AF-{j["ProteinID"]}-F1-model_v4.pdb' # create the filename
        if filename in os.listdir('C:\\Users\\ahmad\\Dropbox\\PC\\Documents\\QMUL\\Research Project\\Proteomics\\Proteomics-extended\\data\\proteome'):   # if the filename is in the human proteome folder
            print(filename)

            df_dssp= dssp_file_df(f'C:\\Users\\ahmad\\Dropbox\\PC\\Documents\\QMUL\\Research Project\\Proteomics\\Proteomics-extended\\data\\proteome\\{filename}\\{filename}') # read in the dssp file
            if j['LastAA_position'] in df_dssp['dssp index']: # if the last AA position is in the dssp file
                relative_asa.append(df_dssp[df_dssp['dssp index'] == j['LastAA_position']]['relative ASA'].values[0]) # append the relative ASA
                secondary_structure.append(df_dssp[df_dssp['dssp index'] == j['LastAA_position']]['secondary structure'].values[0]) # append the secondary structure
                phi.append(df_dssp[df_dssp['dssp index'] == j['LastAA_position']]['phi'].values[0]) # append the phi
                psi.append(df_dssp[df_dssp['dssp index'] == j['LastAA_position']]['psi'].values[0]) # append the psi
                NH01R.append(df_dssp[df_dssp['dssp index'] == j['LastAA_position']]['NH01R'].values[0]) # append the NH01R
                NH01E.append(df_dssp[df_dssp['dssp index'] == j['LastAA_position']]['NH01E'].values[0]) # append the NH01E
                ONH1.append(df_dssp[df_dssp['dssp index'] == j['LastAA_position']]['ONH1'].values[0])   # append the ONH1
                ONHE.append(df_dssp[df_dssp['dssp index'] == j['LastAA_position']]['ONHE'].values[0])   # append the ONHE
                NH2R.append(df_dssp[df_dssp['dssp index'] == j['LastAA_position']]['NH2R'].values[0])   # append the NH2R
                NO2E.append(df_dssp[df_dssp['dssp index'] == j['LastAA_position']]['NO2E'].values[0])   # append the NO2E
                ONH2R.append(df_dssp[df_dssp['dssp index'] == j['LastAA_position']]['ONH2R'].values[0]) # append the ONH2R
                ONH2E.append(df_dssp[df_dssp['dssp index'] == j['LastAA_position']]['ONH2E'].values[0]) # append the ONH2E


            else: # if the last AA position is not in the dssp file
                relative_asa.append('') # append empty string
                secondary_structure.append('')
                phi.append('')
                psi.append('')
                NH01R.append('')
                NH01E.append('')
                ONH1.append('')
                ONHE.append('')
                NH2R.append('')
                NO2E.append('')
                ONH2R.append('')
                ONH2E.append('')
        else:   # if the filename is not in the human proteome folder
            relative_asa.append('') # append empty string
            secondary_structure.append('')  
            phi.append('')
            psi.append('')
            NH01R.append('')
            NH01E.append('')
            ONH1.append('')
            ONHE.append('')
            NH2R.append('')
            NO2E.append('')
            ONH2R.append('')
            ONH2E.append('')


    df['secondary_structure'] = secondary_structure # add the secondary structure to the dataframe
    df['relative ASA'] = relative_asa
    df['phi'] = phi
    df['psi'] = psi
    df['NH01R'] = NH01R
    df['NH01E'] = NH01E
    df['ONH1'] = ONH1
    df['ONHE'] = ONHE
    df['NH2R'] = NH2R
    df['NO2E'] = NO2E
    df['ONH2R'] = ONH2R
    df['ONH2E'] = ONH2E
    return df.to_csv('df_dssp.csv', index=False)   # return the dataframe to a csv file


def residue_depth(file): #This function calculates the residue depth of a given PDB file and makes a dataframe with the results
    parser = PDBParser() #Creates a parser object
    structure = parser.get_structure(f'{file}', file)   #Creates a structure object from the PDB file
    model = structure[0] #Creates a model object from the structure object
    rd = ResidueDepth(model)    #Creates a ResidueDepth object from the model object

    appender = []   #Creates a list to append the results to
    for k in rd.property_keys:  #Iterates through the keys in the ResidueDepth object
        x = rd.property_dict[k] #Creates a variable to store the value of the key
        resdepth = x[0] #Creates a variable to store the value of the first element of the list
        cadepth = x[1]  
        appender.append((resdepth, cadepth))    #Appends the results to the list

    df_rd = pd.DataFrame.from_records(appender, columns=['res_depth', 'ca_depth'])  #Creates a dataframe from the list of results
    df_rd['index'] = range(1, len(df_rd) + 1)   #Adds an index column to the dataframe
    return df_rd




def get_residue_depth_name(df):
    # Function uses proteins in data frame to find all the PDB files in the directory and run the residue depth function on them
    proteinID = []  # List to store the protein IDs
    for p in df['ProteinID']:  # Loop through the proteins in the data frame
        proteinID.append(p)  # Append the protein ID to the list

    proteins = set(proteinID)  # Convert the list to a set to remove duplicates
    print(len(proteins))
    for k in proteins:  # Loop through the set of proteins
        print(k)
        filename = f'AF-{k}-F1-model_v4.pdb'  # Create the filename
        filepath = f'C:\\Users\\ahmad\\Dropbox\\PC\\Documents\\QMUL\\Research Project\\Proteomics\\Proteomics-extended\\data\\rd\\{filename}.csv' #Create the file path
        if os.path.exists(filepath):  # Check if the file already exists in the directory
            print("The file exists!")
        else:
            if filename in os.listdir('C:\\Users\\ahmad\\Dropbox\\PC\\Documents\\QMUL\\Research Project\\Proteomics\\Proteomics-extended\\data\\proteome'):  # Check if the file is in the directory
                df_rd = residue_depth(f'C:\\Users\\ahmad\\Dropbox\\PC\\Documents\\QMUL\\Research Project\\Proteomics\\Proteomics-extended\\data\\proteome\\{filename}\\{filename}')  # Run the residue depth function on the file
                df_rd.to_csv(filepath)  # Save the dataframe to a csv file




def add_values(df): #Function to add the values of the residue depth and cadepth to the dataframe

    res_depth = [] #List to store the residue depth values
    ca_depth = []   
    for i,j in df.iterrows():  #Loop through the dataframe
        file= f'AF-{j["ProteinID"]}-F1-model_v4.pdb.csv'    #Create the filename


        if file in os.listdir('C:\\Users\\ahmad\\Dropbox\\PC\\Documents\\QMUL\\Research Project\\Proteomics\\Proteomics-extended\\data\\rd'):   #Check if the file is in the directory

            df_rd= pd.read_csv(f'C:\\Users\\ahmad\\Dropbox\\PC\\Documents\\QMUL\\Research Project\\Proteomics\\Proteomics-extended\\data\\rd\\{file}')    #Read the csv file
            if j['LastAA_position'] in df_rd['index']:  #Check if the last amino acid position is in the dataframe
                res_depth.append(df_rd.loc[df_rd['index'] == j['LastAA_position']]['res_depth'].values[0])  #Append the residue depth value to the list
                ca_depth.append(df_rd.loc[df_rd['index'] == j['LastAA_position']]['ca_depth'].values[0])    #Append the cadepth value to the list
            else:
                res_depth.append('') #Append a blank value to the list
                ca_depth.append('')  
        else:
            res_depth.append('')    #Append a blank value to the list
            ca_depth.append('') 

    print(res_depth)

    df['res_depth'] = res_depth    #Add the values to the dataframe
    df['ca_depth'] = ca_depth   #Add the values to the dataframe



def get_sasa(file): #file is the pdb file


    parser = PDBParser(QUIET=True) #Creates a parser object
    structure = parser.get_structure(f'{file}', file)   #Creates a structure object from the PDB file

    sasa = ShrakeRupley()
    sasa.compute(structure[0], level="R")

    my_list = [] #Creates a list to store the sasa values
    for chain in structure[0]: #Iterates through the chains in the structure
        for res in chain: #Iterates through the residues in the chain
            my_list.append((res.get_resname(),round(res.sasa,2))) #Appends the residue name and sasa value to the list
 
    aa_df= pd.DataFrame(my_list, columns=['Amino Acids', 'SASA']) #Creates a dataframe from the list
    aa_df['index'] = range(1, len(aa_df) + 1)   #Adds an index column to the dataframe


    amino_acids= {'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', 'ASP':'D', 'ASN':'N', 
    'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y', 'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 
    'MET':'M', 'ALA':'A', 'GLY':'G', 'PRO':'P', 'CYS':'C'} #Creates a dictionary to map the amino acid abbreviation to the AA single letter code

    for key, value in amino_acids.items(): #Iterates through the dictionary
        aa_df['Amino Acids'] = aa_df['Amino Acids'].replace(key, value) #Replaces the abbreviation with the single letter code


    return aa_df #Returns the dataframe


def generating_sasa_columns(df):
    proteinID = [] #List to store the protein IDs
    for p in df['ProteinID']:  #Loop through the proteins in the data frame
        proteinID.append(p) #Append the protein ID to the list
    proteins = set(proteinID)   #Convert the list to a set to remove duplicates
    print(len(proteins))
    for k in proteins:  #Loop through the set of proteins
        
        filename = f'AF-{k}-F1-model_v4.pdb' #Create the filename
        filepath = f'C:\\Users\\ahmad\\Dropbox\\PC\\Documents\\QMUL\\Research Project\\Proteomics\\Proteomics-extended\\data\\sasa\\{filename}.csv'
        if os.path.exists(filepath): #Check if the file already exists in the directory
            print("The file exists!")
        else:
            pdb_filepath = f"C:\\Users\\ahmad\\Dropbox\\PC\\Documents\\QMUL\\Research Project\\Proteomics\\Proteomics-extended\\data\\proteome\\{filename}\\{filename}"
            if filename in os.listdir("C:\\Users\\ahmad\\Dropbox\\PC\\Documents\\QMUL\\Research Project\\Proteomics\\Proteomics-extended\\data\\proteome"): #Check if the file is in the directory
                print(filename)
                aa_df = get_sasa(pdb_filepath)  #Run the get_sasa function on the file
                aa_df.to_csv(filepath)  #Save the dataframe to a csv file
        
    #Code for printing_sasa() function
    sasa_value= [] #List to store the sasa values
    for i,j in df.iterrows():  #Loop through the dataframe
        file= f'AF-{j["ProteinID"]}-F1-model_v4.pdb.csv'    #Create the filename
        
        
        if file in os.listdir("C:\\Users\\ahmad\\Dropbox\\PC\\Documents\\QMUL\\Research Project\\Proteomics\\Proteomics-extended\\data\\sasa"):   #Check if the file is in the directory
            aa_df= pd.read_csv(f'C:\\Users\\ahmad\\Dropbox\\PC\\Documents\\QMUL\\Research Project\\Proteomics\\Proteomics-extended\\data\\sasa\\{file}')  #Read the csv file
            if j['LastAA_position'] in aa_df['index']:  #Check if the last amino acid position is in the dataframe
                sasa_value.append(aa_df.loc[aa_df['index'] == j['LastAA_position']]['SASA'].values[0]) #Append the sasa value to the list
            else:
                sasa_value.append('') #Append a blank value to the list
        else:
            sasa_value.append('')
    print(sasa_value)
    df['SASA']= sasa_value #Add the sasa values to the dataframe
    