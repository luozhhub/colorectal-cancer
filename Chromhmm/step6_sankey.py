import sys, os
import pandas as pd
import plotly.graph_objects as go


def make_windows():
    bedfile_dir = "/home/zhluo/Project/CRC/data_nazhang/step42_chromatine_state_expression/every_week_state"

    for one_file in os.listdir(bedfile_dir):
        if "sorted" not in one_file:
            continue
        else:
            input_bed = os.path.join(bedfile_dir, one_file)

            out_bed = os.path.join("/home/zhluo/Project/CRC/data_nazhang/step43_sankey/segment_200_bin", one_file + ".200bin")
            cmd= "bedtools makewindows -b %s -w 200 -i src >%s" % (input_bed, out_bed)
            os.system(cmd)

def statistic_number():
    bin_file = "/home/zhluo/Project/CRC/data_nazhang/step43_sankey/segment_200_bin/combined.state.binfile.txt"
    df = pd.read_csv(bin_file, sep="\t", names=["chr", "start", "end", "ctrl_state", "week2_state", "week4_state", "week7_state", "week10_state"])
    #print(df[0:5])
    change_state_list = []
    trace_dir = "/home/zhluo/Project/CRC/data_nazhang/step43_sankey/trace/"

    ctrl_week2 = df.groupby(["ctrl_state", "week2_state"]).size().reset_index(name="count")
    ctrl_week2.to_csv(os.path.join(trace_dir, "ctrl_week2.txt"), sep="\t")
    change_state_list.append(ctrl_week2)

    week2_week4 = df.groupby(["week2_state", "week4_state"]).size().reset_index(name="count")
    week2_week4.to_csv(os.path.join(trace_dir, "week2_week4.txt"), sep="\t")
    change_state_list.append(week2_week4)

    week4_week7 = df.groupby(["week4_state", "week7_state"]).size().reset_index(name="count")
    week4_week7.to_csv(os.path.join(trace_dir, "week4_week7.txt"), sep="\t")
    change_state_list.append(week4_week7)

    week7_week10 = df.groupby(["week7_state", "week10_state"]).size().reset_index(name="count")
    week7_week10.to_csv(os.path.join(trace_dir, "week7_week10.txt"), sep="\t")
    change_state_list.append(week7_week10)
    return (change_state_list)


def node_trace():
    #total nodes = 65, 13*5
    lable = ["ctrl_" + "E_" + str(i) for i in range(1,14)] + ["week2_" + "E_" + str(i) for i in range(1,14)] + ["week4_" + "E_" + str(i) for i in range(1,14)] + ["week7_" + "E_" + str(i) for i in range(1,14)]    + ["week10_" + "E_" + str(i) for i in range(1,14)]
    lable_s = ["E" + str(i) for i in range(1,14)] + ["E" + str(i) for i in range(1,14)] + ["E" + str(i) for i in range(1,14)] + ["E" + str(i) for i in range(1,14)]    + ["E" + str(i) for i in range(1,14)]

    print(lable[0:13])
    source,target,value = [], [], []



    ctrl_week2 = pd.read_csv("/home/zhluo/Project/CRC/data_nazhang/step43_sankey/trace/ctrl_week2.txt", sep="\t", index_col=0)
    print(ctrl_week2)
    for i in range(1,14):
        for j in range(1,14):
            if i == j and i != 8:
                continue
            source.append(i)
            target.append(j + 13)
            row = ctrl_week2[(ctrl_week2["ctrl_state"]=="E" + str(i)) & (ctrl_week2["week2_state"]=="E" + str(j))]
            value.append(int(row["count"]))


    week2_week4 = pd.read_csv("/home/zhluo/Project/CRC/data_nazhang/step43_sankey/trace/week2_week4.txt", sep="\t", index_col=0)
    for i in range(1,14):
        for j in range(1,14):
            if i == j and i != 8:
                continue
            source.append(i + 13)
            target.append(j + 13 + 13)
            row = week2_week4[(week2_week4["week2_state"]=="E" + str(i)) & (week2_week4["week4_state"]=="E" + str(j))]
            if row.empty:
                continue
            #print(row["count"])
            value.append(int(row["count"]))


    week4_week7 = pd.read_csv("/home/zhluo/Project/CRC/data_nazhang/step43_sankey/trace/week4_week7.txt", sep="\t", index_col=0)
    for i in range(1,14):
        for j in range(1,14):
            if i == j and i != 8:
                continue
            source.append(i + 13 + 13)
            target.append(j + 13 + 13 + 13)
            row = week4_week7[(week4_week7["week4_state"]=="E" + str(i)) & (week4_week7["week7_state"]=="E" + str(j))]
            if row.empty:
                continue
            value.append(int(row["count"]))


    week7_week10 = pd.read_csv("/home/zhluo/Project/CRC/data_nazhang/step43_sankey/trace/week7_week10.txt", sep="\t", index_col=0)
    for i in range(1,14):
        for j in range(1,14):
            if i == j and i != 8:
                continue
            source.append(i + 13 + 13 + 13 )
            target.append(j + 13 + 13 + 13 + 13)
            row = week7_week10[(week7_week10["week7_state"]=="E" + str(i)) & (week7_week10["week10_state"]=="E" + str(j))]
            if row.empty:
                continue
            value.append(int(row["count"]))

    result = []
    result.append(source)
    result.append(target)
    result.append(value)
    return (result)


def sankey():
    result = node_trace()
    fig = go.Figure(data=[go.Sankey(
    node = dict(
      pad = 15,
      thickness = 20,
      line = dict(color = "black", width = 0.5),
      #label = ["A1", "A2", "B1", "B2", "C1", "C2"],
      color = "blue"
    ),
    link = dict(
      source = result[0], # indices correspond to labels, eg A1, A2, A2, B1, ...
      target = result[1],
      value = result[2]
      ),

    textfont= dict(size=1)
      )])

    fig.update_layout(title_text="Basic Sankey Diagram", font_size=10)
    fig.write_image("./fig1.png", width=600, height=600)








if __name__ == "__main__":
    #make 200 bp windows, the row in each file should be the same
    #make_windows()

    #step43, cd /home/zhluo/Project/CRC/data_nazhang/step43_sankey/segment_200_bin
    #paste -d "\t" /home/zhluo/Project/CRC/data_nazhang/step43_sankey/segment_200_bin/ctrl_13_segments.bed.sorted.200bin /home/zhluo/Project/CRC/data_nazhang/step43_sankey/segment_200_bin/2weeks_13_segments.bed.sorted.200bin /home/zhluo/Project/CRC/data_nazhang/step43_sankey/segment_200_bin/4weeks_13_segments.bed.sorted.200bin /home/zhluo/Project/CRC/data_nazhang/step43_sankey/segment_200_bin/7weeks_13_segments.bed.sorted.200bin /home/zhluo/Project/CRC/data_nazhang/step43_sankey/segment_200_bin/10weeks_13_segments.bed.sorted.200bin | cut -f 1,2,3,4,8,12,16,20 > /home/zhluo/Project/CRC/data_nazhang/step43_sankey/segment_200_bin/combined.state.binfile.txt
    #statistic_number()
    #node_trace()
    sankey()
