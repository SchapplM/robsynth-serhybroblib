% Jacobian time derivative of explicit kinematic constraints of
% palh4m1IC
% with respect to active joint coordinates
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [8x1]
%   Generalized joint coordinates (joint angles)
% qJD [8x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,CB,CE,EP,OT,TA,TD]';
% 
% Output:
% PhiD_a [(no of constraints)x(no. of active joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 23:05
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function PhiD_a = palh4m1IC_kinconstr_impl_act_jacobianD_mdh_sym_varpar(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(8,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [8 1]), ...
  'palh4m1IC_kinconstr_impl_act_jacobianD_mdh_sym_varpar: qJ has to be [8x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [8 1]), ...
  'palh4m1IC_kinconstr_impl_act_jacobianD_mdh_sym_varpar: qJD has to be [8x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'palh4m1IC_kinconstr_impl_act_jacobianD_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_active_jacobianD_matlab.m
% OptimizationMode: 1
% StartTime: 2020-04-11 23:05:40
% EndTime: 2020-04-11 23:05:40
% DurationCPUTime: 0.01s
% Computational Cost: add. (2->2), mult. (6->6), div. (0->0), fcn. (4->4), ass. (0->19)
unknown=NaN(3,5);
t1 = sin(qJ(2));
t3 = sin(qJ(7));
t6 = cos(qJ(2));
t8 = cos(qJ(7));
unknown(1,1) = 0;
unknown(1,2) = -(qJD(2) * t1);
unknown(1,3) = 0;
unknown(1,4) = 0;
unknown(1,5) = -(qJD(7) * t3 * pkin(1));
unknown(2,1) = 0;
unknown(2,2) = (qJD(2) * t6);
unknown(2,3) = 0;
unknown(2,4) = 0;
unknown(2,5) = (qJD(7) * t8 * pkin(1));
unknown(3,1) = 0;
unknown(3,2) = 0;
unknown(3,3) = 0;
unknown(3,4) = 0;
unknown(3,5) = 0;
PhiD_a  = unknown;
