% Jacobian time derivative of explicit kinematic constraints of
% palh1m2IC
% with respect to active joint coordinates
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% qJD [13x1]
%   Generalized joint velocities
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
% 
% Output:
% PhiD_a [(no of constraints)x(no. of active joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:49
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function PhiD_a = palh1m2IC_kinconstr_impl_act_jacobianD_mdh_sym_varpar(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2IC_kinconstr_impl_act_jacobianD_mdh_sym_varpar: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m2IC_kinconstr_impl_act_jacobianD_mdh_sym_varpar: qJD has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2IC_kinconstr_impl_act_jacobianD_mdh_sym_varpar: pkin has to be [20x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_active_jacobianD_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:42:14
% EndTime: 2020-05-02 23:42:14
% DurationCPUTime: 0.02s
% Computational Cost: add. (19->10), mult. (20->15), div. (0->0), fcn. (10->10), ass. (0->9)
t19 = pkin(3) * (qJD(7) + qJD(2));
t18 = pkin(10) * (qJD(3) + qJD(4));
t17 = pkin(1) * qJD(2);
t16 = pkin(5) * qJD(3);
t15 = pkin(6) * qJD(3);
t14 = pkin(17) + qJ(3);
t11 = qJ(3) + qJ(4) + pkin(18);
t10 = qJ(7) + pkin(20) + qJ(2);
t1 = [0, -sin(t10) * t19 - sin(qJ(2)) * t17, 0, 0; 0, cos(t10) * t19 + cos(qJ(2)) * t17, 0, 0; 0, 0, -sin(t14) * t15, 0; 0, 0, cos(t14) * t15, 0; 0, 0, -sin(t11) * t18 - sin(qJ(3)) * t16, 0; 0, 0, cos(t11) * t18 + cos(qJ(3)) * t16, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
PhiD_a = t1;
