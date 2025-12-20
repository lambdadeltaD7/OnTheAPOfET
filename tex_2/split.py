import re
import os

def parse_file(filename):
    with open(filename, 'r', encoding='utf-8') as f:
        content = f.read()
        
    
    # Регулярное выражение для поиска блоков
    pattern = r'([A-Z]+_[A-Z]+)\s*\n(.*?)(?=\n[A-Z]+_[A-Z]+|$)'
    blocks = re.findall(pattern, content, re.DOTALL)
    
    for block_header, block_content in blocks:
    
        # Определяем test_name и dist_name
        test_dist_match = re.match(r'([A-Z]+)_([A-Z]+)', block_header)
        if not test_dist_match:
            continue
            
        test_name, dist_name = test_dist_match.groups()
        
        # Ищем h1 или h2 в блоке
        h_pattern = r'(h[12])\s*=\s*([\d.]+)'
        h_matches = re.findall(h_pattern, block_content)
        
        if not h_matches:
            continue
            
        # Берем первое найденное h значение
        h_type, h_value = h_matches[0]
        
        # Создаем имя файла
        output_filename = f"./notes/res19/splitted/{test_name}_{dist_name}_{h_type}={h_value}.txt"
        
        # Записываем содержимое блока в файл
        with open(output_filename, 'w', encoding='utf-8') as f:
            f.write(f"{block_header}:\n")
            f.write(block_content.strip())


def main():
    input_file = "./notes/res19/all.txt"
    
    if not os.path.exists(input_file):
        print(f"Файл {input_file} не найден!")
        return
    
    parse_file(input_file)


if __name__ == "__main__":
    main()