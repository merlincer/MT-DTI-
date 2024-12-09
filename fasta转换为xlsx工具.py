import tkinter as tk
from tkinter import filedialog, messagebox
from Bio import SeqIO
import pandas as pd
import os

def fasta_to_excel(fasta_file, excel_file):
    """
    将FASTA文件转换为Excel文件。

    参数：
        fasta_file (str): 输入的FASTA文件路径。
        excel_file (str): 输出的Excel文件路径。

    返回：
        None
    """
    # 检查输入文件是否存在
    if not os.path.exists(fasta_file):
        print(f"错误：{fasta_file} 文件不存在。")
        return

    # 读取FASTA文件
    records = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        records.append({"ID": record.id, "Description": record.description, "Sequence": str(record.seq)})

    # 转换为DataFrame
    df = pd.DataFrame(records)

    # 写入Excel文件
    try:
        df.to_excel(excel_file, index=False, sheet_name="FASTA Data")
        print(f"FASTA数据已成功写入 {excel_file}")
    except Exception as e:
        print(f"写入Excel文件时出错：{e}")

def select_input_file():
    """
    打开文件对话框以选择输入的FASTA文件。

    返回：
        str: 选择的文件路径
    """
    file_path = filedialog.askopenfilename(title="选择FASTA文件", filetypes=[("FASTA文件", "*.fasta"), ("所有文件", "*.*")])
    input_entry.delete(0, tk.END)
    input_entry.insert(0, file_path)

def start_conversion():
    """
    开始将FASTA文件转换为Excel文件。
    """
    input_file = input_entry.get()
    if not input_file:
        messagebox.showerror("错误", "请先选择输入文件！")
        return

    output_file = filedialog.asksaveasfilename(title="保存为Excel文件", defaultextension=".xlsx", filetypes=[("Excel文件", "*.xlsx")])
    if not output_file:
        return

    try:
        fasta_to_excel(input_file, output_file)
        messagebox.showinfo("成功", f"文件已成功转换并保存为\n{output_file}")
    except Exception as e:
        messagebox.showerror("错误", f"转换过程中发生错误：{e}")

# 创建主界面
root = tk.Tk()
root.title("FASTA 转 Excel 工具")

# 创建标签和输入框
input_label = tk.Label(root, text="选择FASTA文件：")
input_label.grid(row=0, column=0, padx=10, pady=10)

input_entry = tk.Entry(root, width=50)
input_entry.grid(row=0, column=1, padx=10, pady=10)

browse_button = tk.Button(root, text="浏览", command=select_input_file)
browse_button.grid(row=0, column=2, padx=10, pady=10)

# 创建开始转换按钮
convert_button = tk.Button(root, text="开始转换", command=start_conversion)
convert_button.grid(row=1, column=0, columnspan=3, pady=20)

# 运行主循环
root.mainloop()
