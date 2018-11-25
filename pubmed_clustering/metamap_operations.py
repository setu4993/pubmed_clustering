import json
import os
import string
import subprocess
import tempfile
import logging

translator = str.maketrans(string.punctuation, ' ' * len(string.punctuation))
mm_location = '/opt/public_mm/bin/metamap16'
logger = logging.getLogger()


def documents_disease_list(documents):
    disease_words = []
    document_json_output = []

    for i, document in enumerate(documents):
        json_output = _metamap_fetch(document['all_text'])
        document_json_output.append(json_output)
        disease_words.append(_metamap_matched_words(json_output))
        if i % 10 == 0:
            print(i)

    return disease_words, document_json_output


def _metamap_fetch(sentences, mm_location=mm_location):
    sentences = [sentences]

    if sentences is not None:
        in_file = tempfile.NamedTemporaryFile(mode="wb", delete=False)
        for sentence in sentences:
            in_file.write(b'%r\n' % sentence)
    else:
        raise ValueError("No input defined.")

    in_file.flush()

    out_file = tempfile.NamedTemporaryFile(mode="r", delete=False)

    command = [mm_location]

    command.append("-f")  # number the mappings
    command.append("--negex")  # provide negex output
    command.append("--silent")  # hide header information
    command.append("--prune 25")  # prune the number of candidate concepts constructed
    command.append("--JSONn")  # formatted JSON output
    command.append("-J dsyn,neop")  # restrict to semantic types
    '''
    dsyn - disease or syndrome
    neop - neoplastic process
    '''

    command.append(in_file.name)
    command.append(out_file.name)

    run_metamap = subprocess.Popen(command, stdout=subprocess.PIPE)

    while run_metamap.poll() is None:
        stdout = str(run_metamap.stdout.readline())
        if 'ERROR' in stdout:
            run_metamap.terminate()
            error = stdout.rstrip()
            raise ChildProcessError('Error %s in running MetaMap. Error caught: ' + error)
    output = str(out_file.read())

    os.remove(in_file.name)
    os.remove(out_file.name)

    return output


def _metamap_matched_words(output):
    words_list = []

    try:
        mm_op_json = json.loads(output)
    except ValueError:
        return words_list

    for doc in mm_op_json["AllDocuments"]:
        for utts_num, utts in enumerate(doc["Document"]["Utterances"]):
            for phr_num, phr in enumerate(utts["Phrases"]):
                if phr["Mappings"]:
                    for mappings in phr["Mappings"]:
                        for mapping in mappings["MappingCandidates"]:
                            for word in mapping["CandidatePreferred"].lower().translate(translator).split():
                                if len(word) > 2:
                                    words_list.append(word.lower().translate(translator))
    return words_list
