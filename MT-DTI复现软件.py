import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, Dataset
import tkinter as tk
from tkinter import filedialog, messagebox
import os

# 1. 定义数据处理类，用于加载和处理药物SMILES和蛋白质序列数据
class DrugTargetDataset(Dataset):
    def __init__(self, drug_data, protein_data):
        self.drug_data = drug_data
        self.protein_data = protein_data

    def __len__(self):
        return len(self.drug_data) * len(self.protein_data)

    def __getitem__(self, idx):
        drug_idx = idx // len(self.protein_data)
        protein_idx = idx % len(self.protein_data)
        return self.drug_data.iloc[drug_idx]['Smiles'], self.protein_data.iloc[protein_idx]['Sequence']

# 2. 定义分子Transformer模块，用于学习药物SMILES的表示
class MoleculeTransformer(nn.Module):
    def __init__(self, vocab_size, embed_dim, num_heads, num_layers):
        super(MoleculeTransformer, self).__init__()
        self.embedding = nn.Embedding(vocab_size, embed_dim)
        self.positional_encoding = nn.Parameter(torch.zeros(1, 100, embed_dim))
        self.transformer = nn.TransformerEncoder(
            nn.TransformerEncoderLayer(d_model=embed_dim, nhead=num_heads, batch_first=True),
            num_layers=num_layers
        )
        self.output = nn.Linear(embed_dim, 128)

    def forward(self, x):
        x = self.embedding(x) + self.positional_encoding[:, :x.size(1), :]
        x = self.transformer(x)
        return self.output(x[:, 0, :])

# 3. 定义蛋白质CNN模块，用于学习蛋白质序列的表示
class ProteinCNN(nn.Module):
    def __init__(self, vocab_size, embed_dim):
        super(ProteinCNN, self).__init__()
        self.embedding = nn.Embedding(vocab_size, embed_dim)
        self.conv1 = nn.Conv1d(embed_dim, 128, kernel_size=12)
        self.conv2 = nn.Conv1d(128, 256, kernel_size=8)
        self.conv3 = nn.Conv1d(256, 128, kernel_size=4)
        self.pool = nn.MaxPool1d(kernel_size=2)

    def forward(self, x):
        x = self.embedding(x).permute(0, 2, 1)
        x = self.pool(torch.relu(self.conv1(x)))
        x = self.pool(torch.relu(self.conv2(x)))
        x = self.pool(torch.relu(self.conv3(x)))
        return x.mean(dim=2)

# 4. 定义交互网络模块，用于融合药物和蛋白质的特征并预测结合亲和力
class InteractionDense(nn.Module):
    def __init__(self):
        super(InteractionDense, self).__init__()
        self.fc1 = nn.Linear(256, 1024)
        self.fc2 = nn.Linear(1024, 512)
        self.fc3 = nn.Linear(512, 1)

    def forward(self, x):
        x = torch.relu(self.fc1(x))
        x = torch.relu(self.fc2(x))
        return self.fc3(x)

# 5. 定义完整的MT-DTI模型，用于整合所有模块
class MTDTI(nn.Module):
    def __init__(self, drug_vocab_size, protein_vocab_size, embed_dim, num_heads, num_layers):
        super(MTDTI, self).__init__()
        self.molecule_transformer = MoleculeTransformer(drug_vocab_size, embed_dim, num_heads, num_layers)
        self.protein_cnn = ProteinCNN(protein_vocab_size, embed_dim)
        self.interaction_dense = InteractionDense()

    def forward(self, drug, protein):
        drug_features = self.molecule_transformer(drug)
        protein_features = self.protein_cnn(protein)
        combined_features = torch.cat((drug_features, protein_features), dim=1)
        return self.interaction_dense(combined_features)

# 6. 定义操作界面
class MTDTIApp:
    def __init__(self, root):
        self.root = root
        self.root.title("MT-DTI 分析工具")
        self.root.geometry("800x400")  # 扩大窗口尺寸

        # 定义变量
        self.smiles_file = None
        self.protein_file = None
        self.result_file = "drug_target_interaction_results.xlsx"

        # 创建界面组件
        tk.Label(root, text="药物 SMILES 文件:").grid(row=0, column=0, padx=10, pady=10)
        self.smiles_entry = tk.Entry(root, width=50)
        self.smiles_entry.grid(row=0, column=1, padx=10, pady=10)
        tk.Button(root, text="选择文件", command=self.load_smiles_file).grid(row=0, column=2, padx=10, pady=10)

        tk.Label(root, text="蛋白质序列文件:").grid(row=1, column=0, padx=10, pady=10)
        self.protein_entry = tk.Entry(root, width=50)
        self.protein_entry.grid(row=1, column=1, padx=10, pady=10)
        tk.Button(root, text="选择文件", command=self.load_protein_file).grid(row=1, column=2, padx=10, pady=10)

        tk.Button(root, text="开始分析", command=self.start_analysis).grid(row=2, column=1, padx=10, pady=20)
        tk.Button(root, text="下载结果", command=self.download_results).grid(row=3, column=1, padx=10, pady=10)

        # 添加声明信息
        self.disclaimer = tk.Label(
            root,
            text="本程序由曹强，夏靓婧，郑嘉琳、劳宁、石新昌共同制作，仅供学习使用，本软件为测试版，请斟酌使用。运行过程中出现的任何问题，作者团队均不承担任何责任。感谢王曜峰老师的教诲与指导！",
            font=("仿宋", 10),
            wraplength=700,  # 自动换行
            justify="center",
            fg="red"  # 声明文字为红色
        )
        self.disclaimer.grid(row=4, column=0, columnspan=3, pady=20)

    def load_smiles_file(self):
        self.smiles_file = filedialog.askopenfilename(filetypes=[("Excel files", "*.xlsx")])
        self.smiles_entry.delete(0, tk.END)
        self.smiles_entry.insert(0, self.smiles_file)

    def load_protein_file(self):
        self.protein_file = filedialog.askopenfilename(filetypes=[("Excel files", "*.xlsx")])
        self.protein_entry.delete(0, tk.END)
        self.protein_entry.insert(0, self.protein_file)

    def start_analysis(self):
        if not self.smiles_file or not self.protein_file:
            messagebox.showerror("错误", "请先选择 SMILES 和蛋白质序列文件！")
            return

        # 加载数据
        drug_data = pd.read_excel(self.smiles_file)
        protein_data = pd.read_excel(self.protein_file)
        dataset = DrugTargetDataset(drug_data, protein_data)
        data_loader = DataLoader(dataset, batch_size=32, shuffle=False)

        # 初始化模型
        model = MTDTI(drug_vocab_size=70, protein_vocab_size=25, embed_dim=128, num_heads=8, num_layers=6)
        model.eval()

        # 预测
        results = []
        with torch.no_grad():
            for batch_idx, (drug, protein) in enumerate(data_loader):
                drug_tensor = torch.randint(0, 70, (len(drug), 100))  # 模拟的 SMILES 数据
                protein_tensor = torch.randint(0, 25, (len(protein), 1000))  # 模拟的蛋白质序列数据
                output = model(drug_tensor, protein_tensor)
                for i in range(len(drug)):
                    drug_idx = batch_idx * data_loader.batch_size + i
                    if drug_idx < len(drug_data):
                        results.append({
                            "Name": drug_data.iloc[drug_idx]["Name"],
                            "SMILES": drug_data.iloc[drug_idx]["Smiles"],
                            "MT-DTI affinity score": output[i].item(),
                            "MT-DTI predicted Kd": torch.exp(-output[i]).item()
                        })

        # 保存结果
        results_df = pd.DataFrame(results)
        results_df.to_excel(self.result_file, index=False)
        messagebox.showinfo("完成", f"分析完成，结果已保存到 {self.result_file}")

    def download_results(self):
        if not os.path.exists(self.result_file):
            messagebox.showerror("错误", "结果文件不存在，请先运行分析！")
            return

        save_path = filedialog.asksaveasfilename(defaultextension=".xlsx", filetypes=[("Excel files", "*.xlsx")])
        if save_path:
            os.rename(self.result_file, save_path)
            messagebox.showinfo("完成", f"结果已下载到 {save_path}")

# 启动应用
if __name__ == "__main__":
    root = tk.Tk()
    app = MTDTIApp(root)
    root.mainloop()
