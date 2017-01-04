import pandas as pd
import os

def main():
    with open("../../data/results/simResults.csv","r") as f:
        df = pd.read_csv(f)
    df = df.round(3)
    df.columns = pd.Index([u"Category",u"Mode",u"Network",u"Size",
                           u"Model",u"Treatment Strategy",u"Mean",u"Stddev"],
                          dtype="object")
    networkNames = {"grid":"N1",
                    "scalefree": "N2",
                    "rand": "N3",
                    "crp": "N4"}
    strategyNames = {"none": "No Trt",
                     "myopic": "Myopic",
                     "proximal": "Proximal",
                     "ps": "Policy Search",
                     "all": "Trt All"}

    ## convert values to proper format
    df["Value"] = (df["Mean"].map("{:0.2f}".format)
                   + " (" + df["Stddev"].map("{:0.3f}".format) + ")")
    df["Network"] = df["Network"].map(lambda x : networkNames[x]
                                      if x in networkNames else x)
    df["Treatment Strategy"] = df["Treatment Strategy"].map(
        lambda x : strategyNames[x])

    df = df.drop(["Mean","Stddev"],axis=1)
    grouped = df.groupby(["Category","Mode","Model"])
    for name,group in grouped:
        if name[0] == "toy":
            group = group.drop(["Mode","Model","Category"],axis=1)
            group = pd.pivot_table(group,index=["Network","Size"],
                                   columns="Treatment Strategy",values="Value",
                                   aggfunc=lambda x : " ".join(x))
            group = group.reindex_axis(["No Trt","Proximal","Myopic",
                                        "Trt All","Policy Search"],axis=1)
            # group.columns = group.columns.swaplevel(0,1)
            # group.sortlevel(0,axis=1,inplace=True)
            fileName = os.path.join("../../data/results","_".join(name)+".tex")
            with open(fileName,"w") as f:
                f.write(group.to_latex())
        else:
            group = group.drop(["Mode","Model","Size","Network"],
                               axis=1)
            group = pd.pivot_table(group,index="Category",
                                   columns="Treatment Strategy",values="Value",
                                   aggfunc=lambda x : " ".join(x))
            group = group.reindex_axis(["No Trt","Proximal","Myopic",
                                        "Trt All","Policy Search"],axis=1)
            group.index = group.index.set_names(["WNS"])
            # group.columns = group.columns.swaplevel(0,1)
            # group.sortlevel(0,axis=1,inplace=True)

            fileName = os.path.join("../../data/results","_".join(name)+".tex")
            with open(fileName,"w") as f:
                f.write(group.to_latex(index=False))


if __name__ == "__main__":
    main()
