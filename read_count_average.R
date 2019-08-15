args=commandArgs(T)
read_count_table_path = args[1]
output_dir = args[2]

df_all = read.delim(read_count_table_path, sep = "\t", header=T)
df_avaerage = df_all

df_avaerage[,1] = round(rowMeans(df_all[, 1:3]))
df_avaerage[,2] = round(rowMeans(df_all[, 1:3]))
df_avaerage[,3] = round(rowMeans(df_all[, 1:3]))

df_avaerage[,4] = round(rowMeans(df_all[, 4:6]))
df_avaerage[,5] = round(rowMeans(df_all[, 4:6]))
df_avaerage[,6] = round(rowMeans(df_all[, 4:6]))

df_avaerage[,7] = round(rowMeans(df_all[, 7:9]))
df_avaerage[,8] = round(rowMeans(df_all[, 7:9]))
df_avaerage[,9] = round(rowMeans(df_all[, 7:9]))

df_avaerage[,10] = round(rowMeans(df_all[, 10:12]))
df_avaerage[,11] = round(rowMeans(df_all[, 10:12]))
df_avaerage[,12] = round(rowMeans(df_all[, 10:12]))

df_avaerage[,13] = round(rowMeans(df_all[, 13:15]))
df_avaerage[,14] = round(rowMeans(df_all[, 13:15]))
df_avaerage[,15] = round(rowMeans(df_all[, 13:15]))

if (!file_test("-d", output_dir)){ dir.create(output_dir)}
out_file_path = paste(output_dir, "count_table_average.txt", sep = "/")
write.table(df_avaerage, out_file_path, sep="\t", quote=F, row.names=T, col.names = T, na="NA", eol="\n")