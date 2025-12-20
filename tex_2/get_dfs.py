import pandas as pd
import os

def parse_test_file_simple(filename: str) -> pd.DataFrame:
    """
    Упрощенная версия функции для чтения тестовых файлов
    """
    
    with open(filename, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    # Определяем варьируемый параметр
    h1_count = sum(1 for line in lines if line.strip().startswith('h1 ='))
    h2_count = sum(1 for line in lines if line.strip().startswith('h2 ='))
    
    if h1_count > 1:
        varying_param = 'h1'
        fixed_param = 'h2'
    elif h2_count > 1:
        varying_param = 'h2' 
        fixed_param = 'h1'
    else:
        # Если только одно значение, определяем по имени файла
        varying_param = 'h1' if 'h1=' in filename else 'h2'
        fixed_param = 'h2' if varying_param == 'h1' else 'h1'
    
    # Собираем данные
    data = {}
    current_param_value = None
    n_values = []
    is_et = "ET" in filename
    
    for line in lines:
        line = line.strip()
        
        # Ищем значения варьируемого параметра
        if line.startswith(f"{varying_param} ="):
            current_param_value = float(line.split('=')[1].strip())
            data[current_param_value] = {"emp_powers": [], "asp_power": None}
        
        # Ищем значения n
        elif line.startswith('n ='):
            n_val = line.split('=')[1].split('|')[0].strip()
            if n_val not in n_values:
                n_values.append(n_val)
        
        # Ищем emp_powers
        elif line.startswith('emp_powers ='):
            if current_param_value is not None:
                powers = list(map(float, line.split('=')[1].strip().split()))
                data[current_param_value]["emp_powers"] = powers
        
        # Ищем asp_power (только для ET)
        elif line.startswith('asp_power =') and is_et:
            if current_param_value is not None:
                asp_val = float(line.split('=')[1].strip())
                data[current_param_value]["asp_power"] = asp_val
    
    # Создаем DataFrame
    rows = []
    
    # Добавляем строки для каждого n
    for i, n in enumerate(n_values):
        row = {"index": f"EP (n={n})"}
        for param_val, param_data in sorted(data.items()):
            if i < len(param_data["emp_powers"]):
                row[f"{varying_param}={param_val}"] = param_data["emp_powers"][i]
        rows.append(row)
    
    # Добавляем строку для asp_power если есть
    if is_et:
        asp_row = {"index": "AP"}
        for param_val, param_data in sorted(data.items()):
            if param_data["asp_power"] is not None:
                asp_row[f"{varying_param}={param_val}"] = param_data["asp_power"]
        rows.append(asp_row)
    
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.set_index('index')
    
    return df


for fname in os.listdir("./notes/res19/splitted"):
    df = parse_test_file_simple(f"./notes/res19/splitted/{fname}")
    df.to_csv(f"./notes/res19/csvs/{fname[:fname.index('.txt')]}.csv")