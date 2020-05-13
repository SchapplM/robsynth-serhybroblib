% Jacobian time derivative of explicit kinematic constraints of
% picker2Dm2IC
% with respect to passive joint coordinates
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% qJD [12x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% 
% Output:
% PhiD_p [(no of constraints)x(no. of passive joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 09:21
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function PhiD_p = picker2Dm2IC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2IC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm2IC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar: qJD has to be [12x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2IC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_passive_jacobianD_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 09:20:35
% EndTime: 2020-05-11 09:20:35
% DurationCPUTime: 0.04s
% Computational Cost: add. (42->17), mult. (36->21), div. (0->0), fcn. (18->14), ass. (0->18)
t63 = pkin(2) * (qJD(3) + qJD(9));
t50 = qJD(1) + qJD(2);
t62 = pkin(3) * t50;
t61 = pkin(4) * (qJD(4) + t50);
t60 = pkin(5) * (qJD(1) + qJD(8));
t59 = pkin(1) * qJD(5);
t58 = pkin(6) * qJD(3);
t57 = pkin(6) * qJD(6);
t54 = qJ(1) + qJ(2);
t47 = qJ(4) + t54;
t56 = sin(t47) * t61;
t55 = cos(t47) * t61;
t53 = qJ(1) + qJ(8);
t52 = qJ(3) + qJ(9);
t51 = pkin(8) + qJ(5);
t43 = cos(t52) * t63;
t42 = sin(t52) * t63;
t1 = [-t55 - cos(t54) * t62, 0, -t55, 0, 0, 0, 0, 0, 0, 0; -t56 - sin(t54) * t62, 0, -t56, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, cos(t51) * t59, 0, -cos(t53) * t60, 0, 0, 0, 0; 0, 0, 0, sin(t51) * t59, 0, -sin(t53) * t60, 0, 0, 0, 0; 0, t43 - cos(qJ(3)) * t58, 0, 0, cos(qJ(6)) * t57, 0, t43, 0, 0, 0; 0, t42 - sin(qJ(3)) * t58, 0, 0, sin(qJ(6)) * t57, 0, t42, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
PhiD_p = t1;
