use std::fmt::{Display, Formatter};
use std::vec::IntoIter;

pub struct FastaRecord {
    seqid: Vec<u8>,
    comment: Vec<u8>,
    sequence: Vec<u8>,
}

impl IntoIterator for FastaRecord {
    type Item = u8;
    type IntoIter = IntoIter<u8>;

    fn into_iter(self) -> Self::IntoIter {
        self.sequence.into_iter()
    }
}

impl FastaRecord {
    pub fn new(seqid: &str, comment: &str, buffer_size: usize) -> Self {
        FastaRecord {
            seqid: seqid.as_bytes().to_vec(),
            comment: comment.as_bytes().to_vec(),
            sequence: Vec::with_capacity(buffer_size),
        }
    }

    pub fn push(&mut self, nuc: u8) {
        if nuc != b'\n' {
            self.sequence.push(nuc)
        }
    }
}

impl Display for FastaRecord {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, ">{}", String::from_utf8_lossy(&self.seqid))?;
        if !self.comment.is_empty() {
            writeln!(f, " {}", String::from_utf8_lossy(&self.comment))?;
        } else {
            writeln!(f)?;
        }

        write!(f, "{}", String::from_utf8_lossy(&self.sequence))
    }
}

impl Clone for FastaRecord {
    fn clone(&self) -> Self {
        FastaRecord {
            seqid: self.seqid.clone(),
            comment: self.comment.clone(),
            sequence: self.sequence.clone(),
        }
    }
}


// use crate::contig_naming::sequence_id2str;
// use crate::fasta_records_iterator;
//
// pub struct FastaRecordsMetaIterator {
//     fasta_records_iterator: FastaRecordsMetaIterator,
//     min_genes: usize,
//     buffer_size: usize,
// }
//
//
// impl Iterator for FastaRecordsMetaIterator {
//     type Item = Vec<u8>;
//
//     fn next(&mut self) -> Option<Self::Item> {
//         let mut buffer = Vec::with_capacity(self.buffer_size);
//         let mut prev_record_id = None;
//         let mut found_genes = 0_usize;
//         for record_part in self.fasta_records_iterator {
//             match record_part {
//                 fasta_records_iterator::FastaPartType::ID(id) => {
//                     match prev_record_id.as_ref() {
//                         None => {
//                             let contig = sequence_id2str(id.as_slice());
//                             writeln!(temp_fasta_stream.as_mut().expect("E: Temporary fasta creation failed"), ">{}", contig).expect("E: Failed to write sequence id into temporary fasta file");
//                             prev_record_id = Some(id);
//                             found_genes += 1;
//                         },
//                         Some(prev_id) => {
//                             let replace = if genes {
//                                 if are_genes_of_same_contig(prev_id, &id) {
//                                     // We write it here to avoid duplicating "id" as the compiler fails to prove it will not be moved to a "else" block.
//                                     writeln!(temp_fasta_stream.as_mut().expect("E: Temporary fasta creation failed"), "\n>{}", sequence_id2str(id.as_slice())).expect("E: Failed to write sequence id into temporary fasta file");
//                                     found_genes += 1;
//                                     None // We continue with this gene
//                                 } else {
//                                     Some(id)
//                                 }
//                             } else {
//                                 Some(id)
//                             };
//                             if let Some(replace) = replace {
//                                 if let Some(mut temp_fasta_stream) = temp_fasta_stream.take() {
//                                     temp_fasta_stream.flush().expect("E: Failed to flush temporary fasta file");
//                                     drop(temp_fasta_stream);
//
//                                     if !genes || min_genes == 0 || found_genes >= min_genes {
//                                         let contig = sequence_id2str(if genes {
//                                             get_contig_name(prev_id.as_slice())
//                                         } else {
//                                             prev_id.as_slice()
//                                         });
//                                         let model_name2score = classifier.predict(&seq_file_path, None, buffer_size, &None, false);
//                                         classifier.print_prediction(contig,
//                                                                     &mut output_stream, &model_name2score).unwrap_or_else(|_| panic!("E: Failed to write prediction into '{}", output_file.display()));
//                                     }
//                                 }
//
//                                 temp_fasta_stream = Some(BufWriter::new(File::create(&seq_file_path).expect("E: Failed to create temporary fasta file")));
//                                 {
//                                     found_genes = 1;
//                                     let contig = sequence_id2str(replace.as_slice());
//                                     writeln!(temp_fasta_stream.as_mut().expect("E: Temporary fasta creation failed"), ">{}", contig).expect("E: Failed to write sequence id into temporary fasta file");
//                                 }
//
//                                 let _ = prev_record_id.insert(replace);
//                             }
//                         }
//                     }
//                 },
//                 fasta_records_iterator::FastaPartType::Nuc(nuc) => {
//                     match temp_fasta_stream.as_mut() {
//                         Some(temp_fasta_stream) => write!(temp_fasta_stream, "{}", nuc as char).expect("E: Failed to write into temporary file"),
//                         None => panic!("E: The fasta file '{}' is malformed! a non-description line shows up before description line", fasta.display()),
//                     }
//                 }
//             }
//         }
//     }
// }