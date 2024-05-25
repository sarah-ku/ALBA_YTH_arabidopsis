import urllib
import json
import pandas as pd
import numpy as np
import os
from copy import deepcopy

import tensorflow as tf
from tensorflow.keras.regularizers import l2
from tensorflow.keras.layers import Dense, Input, Add, Activation, ZeroPadding1D, AveragePooling1D, Dropout, Conv2D, Conv1D, MaxPooling1D, LSTM, Flatten, Bidirectional, LayerNormalization, BatchNormalization, GRU
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.models import Model, load_model
from tensorflow.keras import Sequential

from keras.initializers import glorot_uniform
from keras.utils import plot_model
from IPython.display import SVG
import keras.backend as K
from fasta_one_hot_encoder import FastaOneHotEncoder

encoder = FastaOneHotEncoder(
    nucleotides = "ACGT",
    lower = False,
    sparse = False,
    handle_unknown="ignore"
)

def load_data(fold):
  train_one_hot_dir = "model_m6A_model_sequence_fasta_train_SPLIT_" + str(fold) + ".fasta"
  valid_one_hot_dir = "model_m6A_model_sequence_fasta_test_SPLIT_" + str(fold) + ".fasta"

  train_label_dir = "model_m6A_train_label_SPLIT_" + str(fold) + ".json"
  valid_label_dir = "model_m6A_test_label_SPLIT_" + str(fold) + ".json"

  train_one_hot = encoder.transform(path=train_one_hot_dir, verbose=False)
  valid_one_hot = encoder.transform(path=valid_one_hot_dir, verbose=False)

  with open(train_label_dir) as f:
    train_label = json.load(f)
  with open(valid_label_dir) as f:
    valid_label = json.load(f)

  train_one_hot = np.array(train_one_hot)
  valid_one_hot = np.array(valid_one_hot)

  train_label = np.array(train_label)
  valid_label = np.array(valid_label)

  BATCH_SIZE = 512

  train_dataset = (
    tf.data.Dataset
    .from_tensor_slices((train_one_hot,train_label))
    .batch(BATCH_SIZE)
    .shuffle(1024)
  )

  valid_dataset = (
    tf.data.Dataset
    .from_tensor_slices((valid_one_hot,valid_label))
    .batch(BATCH_SIZE)
  )

  return [train_dataset,valid_dataset]


def my_model():
  SEQ_SIZE = 101
  DROP_OUT = 0.25
  L2_PARAM = 0.001
  CONV_WIN = 6

  inputs_OH = Input(shape=(SEQ_SIZE, 4), dtype=tf.float32, name="inputs_OH")

  X = Conv1D(32,kernel_size=CONV_WIN , activation='relu', kernel_regularizer=l2(L2_PARAM), bias_regularizer=l2(L2_PARAM))(inputs_OH)
  X = Dropout(DROP_OUT)(X)
  X = BatchNormalization()(X)

  X = MaxPooling1D(pool_size=2)(X)

  X = Conv1D(32,kernel_size=CONV_WIN , activation='relu', kernel_regularizer=l2(L2_PARAM), bias_regularizer=l2(L2_PARAM))(X)
  X = Dropout(DROP_OUT)(X)
  X = BatchNormalization()(X)

  X = MaxPooling1D(pool_size=2)(X)

  X = Conv1D(64,kernel_size=CONV_WIN , activation='relu', kernel_regularizer=l2(L2_PARAM), bias_regularizer=l2(L2_PARAM))(X)
  X = Dropout(DROP_OUT)(X)
  X = BatchNormalization()(X)

  X = MaxPooling1D(pool_size=2)(X)

  X = Conv1D(64,kernel_size=CONV_WIN , activation='relu', kernel_regularizer=l2(L2_PARAM), bias_regularizer=l2(L2_PARAM))(X)
  X = Dropout(DROP_OUT)(X)
  X = BatchNormalization()(X)

  X = MaxPooling1D(pool_size=2)(X)

 ############################################################

  X = Flatten()(X)

  X = Dense(32, activation='relu', kernel_regularizer=l2(L2_PARAM), bias_regularizer=l2(L2_PARAM))(X)
  X = BatchNormalization()(X)
  X = Dropout(DROP_OUT)(X)

  X = Dense(16, activation='relu', kernel_regularizer=l2(L2_PARAM), bias_regularizer=l2(L2_PARAM))(X)
  X = BatchNormalization()(X)
  X = Dropout(DROP_OUT)(X)

  main_out_1 = tf.keras.layers.Dense(1, activation='sigmoid')(X)

  model = Model(inputs=[inputs_OH], outputs=[main_out_1])

  model.compile(loss=tf.keras.losses.BinaryCrossentropy(), optimizer='adam', metrics=['accuracy'])

  return(model)

def model_fit(model, train_dataset, valid_dataset, N_EPOCH=50):
  model.fit(
      train_dataset,
      validation_data=valid_dataset,
      epochs=N_EPOCH
  )

  yhat = model.predict(valid_dataset)

  return [np.array(yhat),model]

def save_results(yhat,fold):
  model.save('modelseq_MO_fold_' + str(fold) + '.h5')
  np.savetxt(X=yhat,fname='model_m6A_label_predictions_fold_' + str(fold) + '.txt')

for fold in range(5):
  print(fold+1)
  train_dataset, valid_dataset = load_data(fold+1)
  model = my_model()
  print(model.summary())
  yhat, model = model_fit(model, train_dataset, valid_dataset, 50)
  save_results(yhat,fold+1)

weights = model.get_weights()
weights[0].shape
for i in range(40):
  conv_weights = weights[i] 
  print(conv_weights.shape)

print(conv_weights.shape)
  for j in range(conv_weights.shape[2]):
    filter = conv_weights[:,:,j]
    np.savetxt(f"filter_size_64_8_MO_{j}_fold_{fold}.txt", filter)

grads.shape

from keras.models import load_model

for fold in range(5):
  testing_dir = 'model_m6A_model_sequence_fasta_MIC_SLIDE_SPLIT_' + str(fold+1) + '.fasta'
  testing_oh = encoder.transform(path=testing_dir, verbose=False)
  testing_oh = np.array(testing_oh)
  model_path = 'modelseq_MO_fold_' + str(fold+1) + '.h5'
  model = load_model(model_path)
  predictions = model.predict(testing_oh)
  np.savez('PREDICTIONS_m6A_SLIDE_MIC_sequences_' + str(fold + 1) + '.npz', predictions)
