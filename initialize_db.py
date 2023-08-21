import requests, argparse
from tqdm import tqdm
import subprocess as sp


def arguments():
    parser = argparse.ArgumentParser(description="Finds Pantoea plasmids based on homology")
    parser.add_argument('-t', '--threads', required=True, type=str, help="Number of threads to use for building database")
    args = parser.parse_args()
    return args

def download_file_from_google_drive(file_id, destination, chunk_size=1024):
    url = "https://docs.google.com/uc?export=download"

    session = requests.Session()
    params = {'id': file_id, 'confirm': 1}
    response = session.get(url, params=params, stream=True)

    for i, chunk_size_ in save_response_content(response, destination, chunk_size):
        yield i, chunk_size_


def get_confirm_token(response):
    for key, value in response.cookies.items():
        if key.startswith('download_warning'):
            return value

    return None


def save_response_content(response, destination, chunk_size):
    with open(destination, "wb") as f:
        for i, chunk in enumerate(response.iter_content(chunk_size)):
            if chunk:  # filter out keep-alive new chunks
                f.write(chunk)
                yield i, chunk_size



if __name__ == '__main__':
    args=arguments()

    file_id = '19LvSo5rreSe1ycFc0bHhN6-6n3Sc4ort'
    destination = 'db.tar.xz'
    for i, chunk in tqdm(download_file_from_google_drive(file_id, destination), unit='KB', desc='Downloading'):
        pass

    cmd = f"tar -xf {destination}"
    proc1 = sp.run(cmd, shell=True, stdout=sp.PIPE)

    cmd2="makeblastdb -dbtype nucl -in db/chromosome.fasta -out chromdb/chromdb"
    proc2 = sp.run(cmd2, shell=True)

    cmd3="makeblastdb -dbtype nucl -in db/rep_dna.fas -out rep_db/rep_dna"
    proc3 = sp.run(cmd3, shell=True)

    cmd4=f"python build_plas_db.py -i db/original_seq.fasta -o plasmid_db -t {args.threads}"
    proc4 = sp.run(cmd4, shell=True)




