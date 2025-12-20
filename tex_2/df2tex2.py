import os
import pandas as pd
import glob

def parse_csv_files():
    csv_files = glob.glob('./notes/res19/csvs/*.csv')
    
    # Группируем файлы по распределению и фиксированному параметру
    groups = {
        'NORMAL_h1': [],
        'NORMAL_h2': [],
        'CAUCHY_h1': [],
        'CAUCHY_h2': []
    }
    
    # Анализируем имена файлов и читаем данные
    for file_path in csv_files:
        filename = os.path.basename(file_path)
        name_parts = filename.replace('.csv', '').split('_')
        
        if len(name_parts) >= 3:
            test_name = name_parts[0]
            dist_name = name_parts[1]
            fixed_param_full = name_parts[2]  # h1=0 или h2=0
            fixed_param_type = fixed_param_full.split('=')[0]  # h1 или h2
            fixed_param_val = fixed_param_full.split('=')[1]  # 0
            
            key = f"{dist_name}_{fixed_param_type}"
            if key in groups:
                # Читаем CSV файл
                df = pd.read_csv(file_path)
                
                # Добавляем информацию о тесте и фиксированном параметре
                groups[key].append({
                    'test_name': test_name,
                    'fixed_param_val': fixed_param_val,
                    'data': df
                })
    
    return groups

def create_latex_table_for_group(group_data, group_name):
    """Создает LaTeX таблицу для одной группы"""
    
    # Определяем заголовки столбцов из первого файла
    first_df = group_data[0]['data']
    varying_params = [col for col in first_df.columns if col != 'index']
    
    latex_code = f"% Table for {group_name}\n"
    latex_code += "\\begin{table}[!h]\n"
    latex_code += "\\centering\n"
    latex_code += "\\caption{}\n"
    latex_code += "\\begin{tabular}{l c " + " ".join(["c"] * len(varying_params)) + "} \n"
    latex_code += "\\toprule\n"
    latex_code += "test & sample size & " + " & ".join(varying_params) + " \\\\ \n"
    latex_code += "\\midrule\n"
    
    # Для каждого теста в группе
    test_names = sorted(set(item['test_name'] for item in group_data))
    
    for test_name in test_names:
        test_data = [item for item in group_data if item['test_name'] == test_name][0]
        df = test_data['data']
        
        # Добавляем multirow для названия теста
        latex_code += f"\\multirow{{3}}{{*}}{{{test_name}}} \n"
        
        # Для каждой строки размера выборки
        for i, row in df.iterrows():
            size_label = row['index']  # EP (n=100), EP (n=400), EP (n=900)
            if "AP" in size_label:
                continue
            size_label = size_label[ size_label.index('(')+1: size_label.index(')')]
            # Форматируем значения
            values = []
            for param in varying_params:
                value = row[param]
                values.append(f"{value:.2f}".rstrip('0').rstrip('.'))
            
            latex_code += f"& {size_label} & " + " & ".join(values) + " \\\\ \n"
        if test_name!=test_names[-1]:
            latex_code += "\\addlinespace\n\\addlinespace\n"
        else:
            latex_code += "\\bottomrule\n"
    
    latex_code += "\\end{tabular}\n"
    latex_code += "\\end{table}\n"
    return latex_code

def main():
    groups = parse_csv_files()
    
    # Создаем таблицы для каждой группы
    for group_name, group_data in groups.items():
        if group_data:  # Если есть данные для этой группы
            latex_table = create_latex_table_for_group(group_data, group_name)
            with open(f"./notes/res19/tables/{group_name.lower()}_table.tex", 'w', encoding='utf-8') as f:
                f.write(latex_table)
            # print(latex_table)

if __name__ == "__main__":
    main()