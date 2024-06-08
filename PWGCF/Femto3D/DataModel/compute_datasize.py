#!/usr/bin/env python3


from ROOT import gROOT

cleaned_text = []


tables = {}
namespaces = {}

with open("singletrackselector.h") as f:
    for i in f:
        i = i.strip("\n")

        if "namespace" in i and "namespace o2::aod" not in i and "// namespace" not in i:
            # print("Beginning namespace", i)
            namespaces[i] = [i]
            for j in f:
                j = j.strip("\n")
                if j.strip() == "":
                    continue
                namespaces[i].append(j)
                
                if j == namespaces[i][1].replace("{", "}"):
                    break
            for j in namespaces[i]:
                continue
                print(j)
                print("End namespace", namespaces[i])

        if "DECLARE_SOA_TABLE_FULL" in i:
            # print("Beginning table", i)
            full_table = [i]
            for j in f:
                j = j.strip("\n")
                for k in j.strip().split(" "):
                    if k == "":
                        continue
                    full_table.append(k)
                if ")" in j:
                    break
            # print("End table", full_table)
            wait_for_closed = False
            tables[i] = []
            for j in full_table:
                # print(j)
                if "DECLARE" in j:
                    continue
                if "<" in j:
                    if ">" not in j:
                        wait_for_closed = True
                    continue
                if ">" in j:
                    wait_for_closed = False
                    continue
                if wait_for_closed:
                    continue
                tables[i].append(j)
                # print("Variable", f"'{j}'")


for i in namespaces:
    print(i)
    for j in namespaces[i]:
        print(j)
    print("")

if 0:
    for i in tables:
        print(i)
        for j in tables[i]:
            print(j)
        print("")
