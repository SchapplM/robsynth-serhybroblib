% Jacobian time derivative of explicit kinematic constraints of
% palh1m2IC
% with respect to passive joint coordinates
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
% PhiD_p [(no of constraints)x(no. of passive joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:49
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function PhiD_p = palh1m2IC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2IC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m2IC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar: qJD has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2IC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar: pkin has to be [20x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_passive_jacobianD_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:42:14
% EndTime: 2020-05-02 23:42:14
% DurationCPUTime: 0.05s
% Computational Cost: add. (44->20), mult. (36->21), div. (0->0), fcn. (18->14), ass. (0->17)
t52 = pkin(3) * (qJD(7) + qJD(2));
t51 = pkin(8) * (qJD(7) + qJD(10));
t50 = pkin(10) * (qJD(3) + qJD(4));
t49 = pkin(12) * (qJD(8) + qJD(9));
t48 = pkin(2) * qJD(8);
t47 = pkin(4) * qJD(7);
t46 = pkin(7) * qJD(6);
t41 = -qJ(7) + pkin(19);
t32 = -qJ(10) + t41;
t45 = cos(t32) * t51;
t42 = qJ(8) + qJ(9);
t44 = sin(t42) * t49;
t43 = cos(t42) * t49;
t34 = qJ(3) + qJ(4) + pkin(18);
t33 = qJ(7) + pkin(20) + qJ(2);
t30 = sin(t32) * t51;
t1 = [0, -cos(qJ(6)) * t46, -sin(t33) * t52, 0, 0, 0, 0, 0, 0; 0, -sin(qJ(6)) * t46, cos(t33) * t52, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t43 + cos(qJ(8)) * t48, -t43, 0, 0, 0, 0; 0, 0, 0, -t44 + sin(qJ(8)) * t48, -t44, 0, 0, 0, 0; -sin(t34) * t50, 0, -t45 + cos(t41) * t47, 0, 0, -t45, 0, 0, 0; cos(t34) * t50, 0, t30 - sin(t41) * t47, 0, 0, t30, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0;];
PhiD_p = t1;
