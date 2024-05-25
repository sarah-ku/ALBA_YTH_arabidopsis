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

fold = 1
train_one_hot_dir = "model_MO_iclip_model_sequence_fasta_train_fold_SPLIT_" + str(fold) + ".json"
train_one_hot_dir

def load_data_valid(fold):
  valid_one_hot_dir = "model_MO_iclip_model_sequence_fasta_test_fold_SPLIT_" + str(fold) + ".fasta"
  valid_one_hot = encoder.transform(path=valid_one_hot_dir, verbose=False)
  valid_one_hot = np.array(valid_one_hot)
  return valid_one_hot

valid_all = list()
for i in [1,2,3,4,5]:
  valid_all.append(np.array(load_data_valid(i)))

valid_all = np.concatenate(valid_all,axis=0)

valid_all.shape

def load_data(fold):
  train_one_hot_dir = "model_MO_iclip_model_sequence_fasta_train_fold_SPLIT_" + str(fold) + ".fasta"
  valid_one_hot_dir = "model_MO_iclip_model_sequence_fasta_test_fold_SPLIT_" + str(fold) + ".fasta"

  train_label_ALBA_dir = "model_MO_train_label_SPLIT_ALBA_fold_" + str(fold) + ".json"
  valid_label_ALBA_dir = "model_MO_test_label_SPLIT_ALBA_fold_" + str(fold) + ".json"

  train_label_ECT2_dir = "model_MO_train_label_SPLIT_ECT2_fold_" + str(fold) + ".json"
  valid_label_ECT2_dir = "model_MO_test_label_SPLIT_ECT2_fold_" + str(fold) + ".json"

  train_label_both_dir = "model_MO_train_label_SPLIT_both_fold_" + str(fold) + ".json"
  valid_label_both_dir = "model_MO_test_label_SPLIT_both_fold_" + str(fold) + ".json"

  train_one_hot = encoder.transform(path=train_one_hot_dir, verbose=False)
  valid_one_hot = encoder.transform(path=valid_one_hot_dir, verbose=False)

  with open(train_label_ALBA_dir) as f:
    train_label_ALBA = json.load(f)
  with open(valid_label_ALBA_dir) as f:
    valid_label_ALBA = json.load(f)

  with open(train_label_ECT2_dir) as f:
    train_label_ECT2 = json.load(f)
  with open(valid_label_ECT2_dir) as f:
    valid_label_ECT2 = json.load(f)

  with open(train_label_both_dir) as f:
    train_label_both = json.load(f)
  with open(valid_label_both_dir) as f:
    valid_label_both = json.load(f)

  train_one_hot = np.array(train_one_hot)
  valid_one_hot = np.array(valid_one_hot)

  train_label_both = np.array(train_label_both)
  valid_label_both = np.array(valid_label_both)

  train_label_ALBA = np.array(train_label_ALBA)
  valid_label_ALBA = np.array(valid_label_ALBA)

  train_label_ECT2 = np.array(train_label_ECT2)
  valid_label_ECT2 = np.array(valid_label_ECT2)

  BATCH_SIZE = 256

  train_dataset = (
    tf.data.Dataset
    .from_tensor_slices((train_one_hot,(train_label_ECT2,train_label_ALBA,train_label_both)))
    .batch(BATCH_SIZE)
    .shuffle(512)
  )

  valid_dataset = (
    tf.data.Dataset
    .from_tensor_slices((valid_one_hot,(valid_label_ECT2,valid_label_ALBA,valid_label_both)))
    .batch(BATCH_SIZE)
  )
  return [train_dataset,valid_dataset,valid_one_hot]


def my_model():
  SEQ_SIZE = 601
  DROP_OUT = 0.2
  L2_PARAM = 0.001
  CONV_WIN = 6
  FILTER_NUM = 64

  inputs_OH = Input(shape=(SEQ_SIZE, 4), dtype=tf.float32, name="inputs_OH")

  X = Conv1D(64, kernel_size=8,activation=None,kernel_regularizer=l2(L2_PARAM), bias_regularizer=l2(L2_PARAM),name="first_conv")(inputs_OH)
  X = Activation('relu')(X)
  X = Dropout(DROP_OUT)(X)
  X = BatchNormalization()(X)
  X = MaxPooling1D(pool_size=2)(X)

  X = Conv1D(FILTER_NUM,kernel_size=CONV_WIN ,activation=None, kernel_regularizer=l2(L2_PARAM), bias_regularizer=l2(L2_PARAM),name="second_conv")(X)
  X = Activation('relu')(X)
  X = Dropout(DROP_OUT)(X)
  X = BatchNormalization()(X)
  X = MaxPooling1D(pool_size=2)(X)

  X = Conv1D(FILTER_NUM,kernel_size=CONV_WIN , activation='relu', kernel_regularizer=l2(L2_PARAM), bias_regularizer=l2(L2_PARAM))(X)
  X = Dropout(DROP_OUT)(X)
  X = BatchNormalization()(X)
  X = MaxPooling1D(pool_size=2)(X)

  X = Conv1D(FILTER_NUM,kernel_size=CONV_WIN , activation='relu', kernel_regularizer=l2(L2_PARAM), bias_regularizer=l2(L2_PARAM))(X)
  X = Dropout(DROP_OUT)(X)
  X = BatchNormalization()(X)
  X = MaxPooling1D(pool_size=2)(X)

  X = Conv1D(FILTER_NUM,kernel_size=CONV_WIN , activation='relu', kernel_regularizer=l2(L2_PARAM), bias_regularizer=l2(L2_PARAM))(X)
  X = Dropout(DROP_OUT)(X)
  X = BatchNormalization()(X)
  X = MaxPooling1D(pool_size=2)(X)

  ############################################################

  X = Flatten()(X)

  O1 = Dense(32, activation='relu', kernel_regularizer=l2(L2_PARAM), bias_regularizer=l2(L2_PARAM))(X)
  O1 = BatchNormalization()(O1)
  O1 = Dropout(DROP_OUT)(O1)

  O2 = Dense(32, activation='relu', kernel_regularizer=l2(L2_PARAM), bias_regularizer=l2(L2_PARAM))(X)
  O2 = BatchNormalization()(O2)
  O2 = Dropout(DROP_OUT)(O2)

  O3 = Dense(32, activation='relu', kernel_regularizer=l2(L2_PARAM), bias_regularizer=l2(L2_PARAM))(X)
  O3 = BatchNormalization()(O3)
  O3 = Dropout(DROP_OUT)(O3)

  main_out_1 = tf.keras.layers.Dense(1, activation='sigmoid')(O1)
  main_out_2 = tf.keras.layers.Dense(1, activation='sigmoid')(O2)
  main_out_3 = tf.keras.layers.Dense(1, activation='sigmoid')(O3)

  model = Model(inputs=[inputs_OH], outputs=[main_out_1,main_out_2,main_out_3])

  model.compile(loss=tf.keras.losses.BinaryCrossentropy(), optimizer='adam', metrics=['accuracy'])

  return(model)



def model_fit(model, train_dataset, valid_dataset, N_EPOCH=50):
  model.fit(
      train_dataset,
      validation_data=valid_dataset,
      epochs=N_EPOCH
  )

  yhat = model.predict(valid_all)

  ta = tf.cast(valid_one_hot,tf.float32)

  inp_tensor_list = [ta]

  with tf.GradientTape() as tape:
      tape.watch(ta)
      pred = model(inp_tensor_list)

  grads_1 = tape.gradient(pred[0], inp_tensor_list)

  with tf.GradientTape() as tape:
      tape.watch(ta)
      pred = model(inp_tensor_list)

  grads_2 = tape.gradient(pred[1], inp_tensor_list)

  with tf.GradientTape() as tape:
      tape.watch(ta)
      pred = model(inp_tensor_list)

  grads_3 = tape.gradient(pred[2], inp_tensor_list)

  grads_1 = np.array(grads_1)
  grads_2 = np.array(grads_2)
  grads_3 = np.array(grads_3)

  return [np.array(yhat),grads_1,grads_2,grads_3,model]

def save_results(yhat,model,features1,grads_1,grads_2,grads_3,fold):
  model.save('modelseq_MO_5CV_fold_' + str(fold) + '.h5')

  np.savez("PREDICTIONS_MO_5CV_fold_" + str(fold) + ".npz", yhat)

  np.savez("GRADS_1_MO_5CV_fold_" + str(fold) + ".npz", grads_1)
  np.savez("GRADS_2_MO_5CV_fold_" + str(fold) + ".npz", grads_2)
  np.savez("GRADS_3_MO_5CV_fold_" + str(fold) + ".npz", grads_3)

  np.savez("FEATURES_1_5CV_fold_" + str(fold) + ".npz", features1)

  weights = model.get_weights()

  conv_weights = weights[0]  
  print(conv_weights.shape)
  np.savez("WEIGHTS_CONV_1_5CV_fold_" + str(fold) + ".npz", conv_weights)


for fold in range(5):
  print(fold+1)
  train_dataset, valid_dataset, valid_one_hot = load_data(fold+1)
  model = my_model()
  print(model.summary())
  yhat, grads_1,grads_2,grads_3, model = model_fit(model, train_dataset, valid_dataset, 18)

  print("Feature extraction from the model")
  feature_extractor = tf.keras.Model(
   inputs=model.inputs,
   outputs=[model.get_layer(name="first_conv").output],
  )

  features = feature_extractor(valid_all)
  print(np.array(features).shape)
  features_layer1 = features
  #features_layer2 = features[1]
  save_results(yhat,model,np.array(features_layer1),grads_1[0],grads_2[0],grads_3[0],fold+1)
