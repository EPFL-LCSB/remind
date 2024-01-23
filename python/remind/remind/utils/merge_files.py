import pandas as pd
import glob
import time

'to merge files e.g for different growth rates etc and concatanete them'





def merge_files(allFiles):
    list_ = []
    frame = pd.DataFrame()
    for file_ in allFiles:
        df = pd.read_hdf(file_)
        list_.append(df)
    frame = pd.concat(list_, ignore_index=True)
    return frame



'input whatever extension you have-- if csv put ".csv" etc ....'
'read in chunks in case you ll have too many files'
def merge_files_w_chunks(path,extension=".h5",n=100):
    allFiles = glob.glob(path + "/*"+extension)
    #n=100
    chunks=[allFiles[i:i + n] for i in range(0, len(allFiles), n)]

    list_=[]
    for i in range(len(chunks)):
        allFiles=chunks[i]
        print(i)
        df=merge_files(allFiles)
        list_.append(df)

    frame = pd.concat(list_, ignore_index=True)

    return frame


def merge_files_for_sliced_yield(path,model,extension=".h5"):
    allFiles = glob.glob(path + "/*"+extension)
    list_ = []
    frame = pd.DataFrame()
    for file_ in allFiles:
        if extension==".h5":
            df = pd.read_hdf(file_)
        if extension==".csv":
            df = pd.read_csv(file_)
        #get_subs_uptake(df,model)
        list_.append(df)
    frame = pd.concat(list_, ignore_index=True)
    return frame


