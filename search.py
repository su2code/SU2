import os

def search_keyword_in_file(file_path, keyword):
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            for line in file:
                if keyword in line:
                    return True
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
    return False

def search_keyword_in_directory(directory, keyword):
    matches = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(('.cpp', '.hpp', '.inl')):
                file_path = os.path.join(root,file)
                if search_keyword_in_file(file_path,keyword):
                    matches.append(file_path)
    
    return matches

if __name__ == "__main__":
    keyword_to_search = input("Enter the keyword: ")
    directory_to_search = os.getcwd()
    results = search_keyword_in_directory(directory_to_search,keyword_to_search)
    if results:
        print(f"'{keyword_to_search}' found in the following files: ") 
        for file_path in results:
            print(f" - {file_path}")
    else:
        print(f"No matches found for '{keyword_to_search}'.")
