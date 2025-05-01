## this is to create a meta data for the 493 voxels in the MOp dataset

voxel.type <- rep('Tissue',493)
empty_index <- c(8:17, 30:34)
voxel.type[empty_index] <- "Empty"
table(voxel.type)

Voxel <- paste0("R", 1:493)

Row <- rep(1:17, 29)
Column <- rep(1:29, each = 17)


x <- data.frame(Voxel = Voxel,
                Row = Row,
                Column = Column,
                voxel.type = voxel.type)
write.csv(x,"MOp_meta_data.csv")
