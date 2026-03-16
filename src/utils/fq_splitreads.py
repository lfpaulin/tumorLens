import sys


# TODO: add +- 1.2kb on each side
class FastqSAM:
    def __init__(self, new_read_len=1000):
        # expects 3 tab-sep cols: read_name, seq, qual
        # this is done by samtools view FASTQ-file | cut -f 1,10,11
        self.new_read_len = new_read_len
        self.read_name, self.seq, self.seq_quality, self.read_len = "", "", 0, 0

    def process_read(self, line):
        self.set_read_info(line)
        self.split_read()

    def set_read_info(self, line):
        [self.read_name, self.seq, self.seq_quality] = line.split("\t")
        self.read_len = len(self.seq)

    def split_read(self):
        for start_pos in range(0, self.read_len, self.new_read_len):
            end_pos = start_pos+self.new_read_len if start_pos+self.new_read_len < self.read_len else self.read_len
            self.make_fastq(start_pos, end_pos)

    def make_fastq(self, read_start, read_end):
        new_read = self.seq[read_start:read_end]
        new_qual = self.seq_quality[read_start:read_end]
        if len(new_read) != 0:
            print(f'@{self.read_name}_{read_start+1}_{read_end}\n{new_read}\n+\n{new_qual}', end="\n")


def main():
    # read_length_target = 1000  # 1kb
    read_length_target = 3000  # 3kb
    fastq_read = FastqSAM(read_length_target)
    for line in sys.stdin:
        fastq_read.process_read(line.rstrip("\n"))


# main
if __name__ == '__main__':
    main()
