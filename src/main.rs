// split the peptide in k-mers with a window size of 1
fn split_sequence(seq: &str, k: usize) -> Vec<String> {
    let mut kmers = Vec::new();
    let mut i = 0;
    while i < seq.len() - k + 1 {
        kmers.push(seq[i..i + k].to_string());
        i += 1;
    }
    kmers
}

fn main() {
    let peptide = "YLLDLHSYL";
    let kmers = split_sequence(peptide, 3);
    println!("{:?}", kmers);    
}