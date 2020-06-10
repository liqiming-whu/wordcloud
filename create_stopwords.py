#!/usr/bin/env python3
import os
import re
from nltk.tokenize import word_tokenize
from nltk.probability import FreqDist
from get_abstract import get_abs

def clean_str(text):
    text = re.sub(r"[^A-Za-z]", " ", text)
    text = re.sub(r"\s{2,}", " ", text)

    return text


def read_abstract(path):
    text = ""
    for file in os.listdir(path):
        print(file)
        text += get_abs(os.path.join(path, file))
    return text


def get_stopwords(text):
    words = word_tokenize(clean_str(text))
    fdist = FreqDist(words)
    stopwords = fdist.most_common(200)
    print(stopwords)
    with open("stopwords.txt", "w", encoding="utf-8") as f:
        for word in stopwords:
            print(word[0], file=f)


def main():
    get_stopwords(read_abstract("download"))


if __name__ == "__main__":
    main()