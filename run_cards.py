"""
    This program uses card_runner.exe to run all the cards in the ./cards directory

    Matias Mantiñan
"""

import glob
import subprocess
import os

def print_progress_bar(iteration, total, prefix='', suffix='', decimals=1, length=50, fill='█'):
    """
    Call in a loop to create a progress bar.
    @param iteration: Required: current iteration (int)
    @param total: Required: total iterations (int)
    @param prefix: Optional: prefix string (str)
    @param suffix: Optional: suffix string (str)
    @param decimals: Optional: positive number of decimals in percent complete (int)
    @param length: Optional: character length of bar (int)
    @param fill: Optional: bar fill character (str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end='\r')
    # Print New Line on Complete
    if iteration == total: 
        print()



def main():
    card_runner = "./bin/card_runner.exe"
    path = "./cards"
    output_path = "./output"
    # Create the directory if it doesn't exist
    os.makedirs(output_path, exist_ok=True)

    nevents = "10000"

    cards_list = []
    for card in glob.glob(path+"/*.cmnd"):
        cards_list.append(card)


    ignore_file = "./cards/IGNORE"
    ignore_cards = []

    with open(ignore_file) as file:
        for line in file:
            ignore_cards.append(line.rstrip())

    num_cards = len(cards_list)

    for card_index,card in enumerate(cards_list):
        # skip cards that are in the ignore file
        print(card)
        if card.split("/")[-1] in ignore_cards:
            continue
        
        print_progress_bar(card_index + 1, num_cards, prefix='Progress:', suffix='Complete', length=50)
        
        output_file_name = card.split("/")[-1] # get only the card name
        output_file_name = output_file_name[:-5] # erase the .cmnd extension, last 5 characters
        output_file_name = output_path +"/" +output_file_name

        command = [card_runner,card,output_file_name,nevents]

        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # The process is now running in the background
        #print(f"Started background process with PID: {process.pid}")

        process.wait()






if __name__ == "__main__":
    main()
