use bio::io::fasta;
use clap::{App, Arg};
use rusqlite;


// read in in FASTA file and return a vector of sequences 
fn read_fasta(filename: &str) -> Vec<(String, usize)> {
    let mut i: usize = 1;
    let mut seqs = Vec::new();
    let reader = fasta::Reader::from_file(filename).unwrap();
    for result in reader.records() {
        let record = result.unwrap();
        let seq_str = std::str::from_utf8(record.seq()).unwrap();
        seqs.push((seq_str.to_string(), i));
        i += 1;
    }
    seqs
}

// split the peptide into k-mers with a window size of 1 and store also the index of that k-mer
fn split_sequence(seq: &str, k: usize) -> Vec<(String, usize)> {
    let mut kmers = Vec::new();
    let mut i: usize = 0;
    while i < seq.len() - k + 1 {
        kmers.push((seq[i..i + k].to_string(), i));
        i += 1;
    }
    kmers
}

// connect to SQLite DB, call it proteome.db
fn connect() -> rusqlite::Connection {
    rusqlite::Connection::open("proteome.db").unwrap()
}

// create a kmers --> index table in the DB
fn create_kmers_table(conn: &rusqlite::Connection) {
    conn.execute(
        "CREATE TABLE IF NOT EXISTS kmers (
                  kmer            TEXT NOT NULL,
                  idx             INTEGER NOT NULL
                  )",
        rusqlite::params![],
    )
    .unwrap();
}

// insert kmers into the table
fn insert_kmers(conn: &mut rusqlite::Connection, kmers: &[(String, usize)], protein_count: &usize) {
    let tx = conn.transaction().unwrap();
    let mut stmt = tx
        .prepare("INSERT INTO kmers (kmer, idx) VALUES (?1, ?2)")
        .unwrap();

    for kmer in kmers {
        stmt.execute(rusqlite::params![kmer.0, (protein_count * 1000000) + kmer.1])
            .unwrap();
    }

    drop(stmt); // Explicitly drop stmt before committing the transaction

    // create index on kmer column
    tx.execute("CREATE INDEX IF NOT EXISTS kmer_idx ON kmers (kmer)", rusqlite::params![])
        .unwrap();

    tx.commit().unwrap();
}


fn main() {
    let matches = App::new("Preprocess proteome.")
        .arg(
            Arg::with_name("proteome")
                .short('p')
                .long("proteome")
                .value_name("FILE")
                .help("Input FASTA file")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("k")
                .short('k')
                .long("k_value")
                .value_name("K")
                .help("Value of k for k-mers")
                .takes_value(true)
                .required(true),
        )
        .get_matches();

    let filename = matches.value_of("proteome").unwrap();
    let k: usize = matches
        .value_of("k")
        .unwrap()
        .parse()
        .unwrap_or_else(|_| {
            eprintln!("Error: k must be an integer");
            std::process::exit(1);
        });

    let seqs = read_fasta(filename);
    let mut conn = connect();
    create_kmers_table(&conn);

    for seq in seqs {
        let kmers = split_sequence(&seq.0, k);
        insert_kmers(&mut conn, &kmers, &seq.1);
    }
}