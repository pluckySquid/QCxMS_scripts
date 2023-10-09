import argparse

def main():
    parser.add_argument('csv_file', type=str, help='csv_file')

    args = parser.parse_args()

    csv_file = args.csv_file

    

if __name__ == "__main__":
    main()